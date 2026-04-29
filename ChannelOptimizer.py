from itertools import combinations

import numpy as np

from ComponentDatabase import evaluate_filter, evaluate_sensor, evaluate_source


def material_names_from_stack(material_library, web_layers=None):
    if web_layers:
        names = []
        for layer in web_layers:
            name = layer["mat_var"].get()
            if name != "Air" and name not in names:
                names.append(name)
        if names:
            return names
    return [name for name in material_library.keys() if name != "Air"]


def build_effective_alpha_matrix(wavelengths, material_library, materials, source_spectrum, channels):
    matrix = np.zeros((len(channels), len(materials)), dtype=float)
    weights = np.zeros(len(channels), dtype=float)

    for row_idx, channel in enumerate(channels):
        weight = source_spectrum * channel["filter_spectrum"] * channel["sensor_spectrum"]
        weight_area = np.trapezoid(weight, wavelengths)
        weights[row_idx] = weight_area

        if weight_area <= 0:
            continue

        for col_idx, material in enumerate(materials):
            alpha = np.asarray(material_library[material]["alpha"], dtype=float)
            matrix[row_idx, col_idx] = np.trapezoid(weight * alpha, wavelengths) / weight_area

    return matrix, weights


def score_matrix(matrix, channel_weights):
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    nonzero = singular_values[singular_values > 1e-12]
    rank = int(np.linalg.matrix_rank(matrix, tol=1e-9))

    if len(nonzero) >= 2:
        condition = float(nonzero[0] / nonzero[-1])
        conditioning_score = float(nonzero[-1] / nonzero[0])
    elif len(nonzero) == 1:
        condition = float("inf")
        conditioning_score = 0.0
    else:
        condition = float("inf")
        conditioning_score = 0.0

    column_norms = np.linalg.norm(matrix, axis=0)
    usable_columns = column_norms > 1e-12
    if np.count_nonzero(usable_columns) >= 2:
        normalized = matrix[:, usable_columns] / column_norms[usable_columns]
        corr = np.abs(normalized.T @ normalized)
        off_diag = corr[~np.eye(corr.shape[0], dtype=bool)]
        max_column_similarity = float(np.max(off_diag)) if off_diag.size else 0.0
    else:
        max_column_similarity = 1.0

    orthogonality = max(0.0, 1.0 - max_column_similarity)
    rank_fraction = rank / matrix.shape[1] if matrix.shape[1] else 0.0
    active_weight_fraction = float(np.mean(channel_weights > 1e-9)) if len(channel_weights) else 0.0

    score = (
        1000.0 * rank_fraction
        + 300.0 * orthogonality
        + 200.0 * conditioning_score
        + 50.0 * active_weight_fraction
    )

    return {
        "score": score,
        "rank": rank,
        "condition": condition,
        "orthogonality": orthogonality,
        "max_column_similarity": max_column_similarity,
        "singular_values": singular_values,
    }


def _rank_partial_sets(channel_count, score_partial):
    return channel_count >= 1 and score_partial


def rank_orthogonal_combinations(
    wavelengths,
    material_library,
    component_database,
    materials=None,
    channel_count=None,
    top_n=20,
    beam_width=250,
    min_channel_weight=1e-6,
):
    if materials is None:
        materials = [name for name in material_library.keys() if name != "Air"]
    if channel_count is None:
        channel_count = len(materials)

    ranked_results = []
    sources = component_database.get("sources", [])
    filters = component_database.get("filters", [])
    sensors = component_database.get("sensors", [])

    for source_def in sources:
        source_spectrum = evaluate_source(source_def, wavelengths)
        candidate_channels = []

        for filter_def in filters:
            filter_spectrum = evaluate_filter(filter_def, wavelengths)
            for sensor_def in sensors:
                sensor_spectrum = evaluate_sensor(sensor_def, wavelengths)
                weight_area = np.trapezoid(source_spectrum * filter_spectrum * sensor_spectrum, wavelengths)
                if weight_area < min_channel_weight:
                    continue
                candidate_channels.append({
                    "filter": filter_def,
                    "sensor": sensor_def,
                    "filter_spectrum": filter_spectrum,
                    "sensor_spectrum": sensor_spectrum,
                    "weight_area": weight_area,
                    "filter_key": filter_def.get("name", f"{filter_def.get('center_nm')}_{filter_def.get('fwhm_nm')}"),
                })

        if len(candidate_channels) < channel_count:
            continue

        beam = [()]
        for depth in range(channel_count):
            expanded = []
            for partial in beam:
                start = partial[-1] + 1 if partial else 0
                for next_idx in range(start, len(candidate_channels)):
                    next_filter_key = candidate_channels[next_idx]["filter_key"]
                    if any(candidate_channels[idx]["filter_key"] == next_filter_key for idx in partial):
                        continue
                    combo = partial + (next_idx,)
                    channels = [candidate_channels[idx] for idx in combo]
                    matrix, weights = build_effective_alpha_matrix(
                        wavelengths,
                        material_library,
                        materials,
                        source_spectrum,
                        channels,
                    )
                    metrics = score_matrix(matrix, weights)
                    if _rank_partial_sets(depth + 1, True):
                        expanded.append((metrics["score"], combo))

            expanded.sort(reverse=True, key=lambda item: item[0])
            beam = [combo for _, combo in expanded[:beam_width]]

        for combo in beam:
            channels = [candidate_channels[idx] for idx in combo]
            matrix, weights = build_effective_alpha_matrix(
                wavelengths,
                material_library,
                materials,
                source_spectrum,
                channels,
            )
            metrics = score_matrix(matrix, weights)
            ranked_results.append({
                "source": source_def,
                "channels": channels,
                "materials": materials,
                "matrix": matrix,
                "weights": weights,
                **metrics,
            })

    ranked_results.sort(key=lambda item: item["score"], reverse=True)
    return ranked_results[:top_n]


def rank_exhaustive_combinations(
    wavelengths,
    material_library,
    component_database,
    materials,
    channel_count,
    top_n=20,
    min_channel_weight=1e-6,
):
    results = []
    for source_def in component_database.get("sources", []):
        source_spectrum = evaluate_source(source_def, wavelengths)
        candidate_channels = []

        for filter_def in component_database.get("filters", []):
            filter_spectrum = evaluate_filter(filter_def, wavelengths)
            for sensor_def in component_database.get("sensors", []):
                sensor_spectrum = evaluate_sensor(sensor_def, wavelengths)
                weight_area = np.trapezoid(source_spectrum * filter_spectrum * sensor_spectrum, wavelengths)
                if weight_area >= min_channel_weight:
                    candidate_channels.append({
                        "filter": filter_def,
                        "sensor": sensor_def,
                        "filter_spectrum": filter_spectrum,
                        "sensor_spectrum": sensor_spectrum,
                        "weight_area": weight_area,
                        "filter_key": filter_def.get("name", f"{filter_def.get('center_nm')}_{filter_def.get('fwhm_nm')}"),
                    })

        for combo in combinations(candidate_channels, channel_count):
            filter_keys = [channel["filter_key"] for channel in combo]
            if len(set(filter_keys)) != len(filter_keys):
                continue
            matrix, weights = build_effective_alpha_matrix(
                wavelengths,
                material_library,
                materials,
                source_spectrum,
                list(combo),
            )
            metrics = score_matrix(matrix, weights)
            results.append({
                "source": source_def,
                "channels": list(combo),
                "materials": materials,
                "matrix": matrix,
                "weights": weights,
                **metrics,
            })

    results.sort(key=lambda item: item["score"], reverse=True)
    return results[:top_n]
