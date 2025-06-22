// initialize temperature profile
initialize(temperature_previous);

while (worst_dt > temp_tolerance && iteration <= max_iterations) {
    // finite difference
    FOR_ALL(i, 1, height + 1,