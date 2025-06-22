// initialize temperature profile
initialize(temperature_previous);

while (worst_dt > temp_tolerance && iteration <= max_iterations) {
    // finite difference
    for (i = 1; i < height + 1; i++) {
        for (j = 1; j < width + 1; j++) {
            temperature(i, j) = 0.25 * (  temperature_previous(i + 1, j)
                                    + temperature_previous(i - 1, j)
                                    + temperature_previous(i, j + 1)
                                    + temperature_previous(i, j - 1));
        }
    }
}