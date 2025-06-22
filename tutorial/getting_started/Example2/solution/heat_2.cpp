// initialize temperature profile
initialize(temperature_previous);

while (worst_dt > temp_tolerance && iteration <= max_iterations) {
    // finite difference
    FOR_ALL(i, 1, height + 1,
            j, 1, width + 1, {
        temperature(i, j) = 0.25 * (  temperature_previous(i + 1, j)
                                    + temperature_previous(i - 1, j)
                                    + temperature_previous(i, j + 1)
                                    + temperature_previous(i, j - 1));
    });
}