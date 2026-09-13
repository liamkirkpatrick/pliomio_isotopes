%% Run a provisional end-to-end SWIM test with Allan Hills isotope data
%
% This script:
%   1. reads tests/test_data.csv;
%   2. generates a SWIM forward state space with the candidate legacy baseline;
%   3. reconstructs condensation, source, and surface temperatures;
%   4. exports one detailed T0 = 10 C to Tcond = -30 C trajectory;
%   5. writes Python-readable numeric fixtures and reproducibility metadata
%      to a timestamped directory under tests/generated/.
%
% Run from any working directory with:
%   run('tests/run_swim_allan_hills_test.m')
%
% This is a provisional integration fixture, not yet a frozen parity baseline.

clearvars;

%% Paths

% MATLAB returns mfilename('fullpath') without the .m extension for scripts.
script_path = [mfilename('fullpath') '.m'];
tests_dir = fileparts(script_path);
project_root = fileparts(tests_dir);
legacy_dir = fullfile(project_root, 'legacy_matlab');
input_path = fullfile(tests_dir, 'test_data.csv');

assert(isfolder(legacy_dir), 'Missing legacy MATLAB directory: %s', legacy_dir);
assert(isfile(input_path), 'Missing Allan Hills input file: %s', input_path);

addpath(legacy_dir);

%% Test configuration

% A 1 deg C grid keeps this end-to-end test reasonably quick. Change both
% increments to 0.1 deg C to reproduce the legacy state-space resolution.
config.Tsite = [-70, 0, 1];
config.Tsource = [0, 28, 1];
config.RHsource = [];
config.supersaturation_a = 1;
config.supersaturation_b = 0.00525;
config.supersaturation_c = 0;
config.closure = 'local';
config.reanalysis = 'ncep';
config.hemisphere = 1;  % Southern Hemisphere
config.season = 'annual';
config.reconstruction_method = 1;  % log-d18O plus dln
config.internal_trajectory_step_degC = 0.1;  % hard-coded by legacy wrapper
config.trajectory.Tsource_C = 10;
config.trajectory.Tcond_C = -30;
config.trajectory.step_degC = 0.1;

started_at = datetime('now', 'TimeZone', 'UTC');
started_at_for_id = started_at;
started_at_for_id.Format = 'yyyyMMdd''T''HHmmss''Z''';
run_id = ['allan_hills_swim_' char(started_at_for_id)];
output_dir = fullfile(tests_dir, 'generated', run_id);
mkdir(output_dir);

fprintf('SWIM Allan Hills test run: %s\n', run_id);
fprintf('Output directory: %s\n', output_dir);

%% Read and validate observations

input_data = readtable(input_path, 'VariableNamingRule', 'preserve');
required_variables = {'sample_id', 'd18O_smow', 'dD_smow'};
missing_variables = setdiff(required_variables, input_data.Properties.VariableNames);
assert(isempty(missing_variables), ...
    'Input CSV is missing required columns: %s', strjoin(missing_variables, ', '));

observed_d18O = input_data.d18O_smow;
observed_dD = input_data.dD_smow;
valid_rows = isfinite(observed_d18O) & isfinite(observed_dD) & ...
    (1 + observed_d18O ./ 1000 > 0) & (1 + observed_dD ./ 1000 > 0);
assert(any(valid_rows), 'Input CSV contains no valid d18O/dD observation pairs.');

%% Capture provenance before running the model

metadata.schema_version = 1;
metadata.test_name = 'Allan Hills SWIM end-to-end test';
metadata.baseline_status = 'provisional';
metadata.run_id = run_id;

started_at_for_text = started_at;
started_at_for_text.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX';
metadata.started_at_utc = char(started_at_for_text);

metadata.project.git_commit = git_value(project_root, 'rev-parse HEAD');
git_status = git_value(project_root, 'status --porcelain');
metadata.project.git_dirty = ~isempty(git_status);
metadata.project.git_status = git_status;
metadata.project.legacy_upstream_commit = ...
    '4db23b79aab8f2154254d9111f77e70a5ce8427d';

metadata.matlab.version = version;
metadata.matlab.release = version('-release');
metadata.matlab.computer = computer;
metadata.matlab.architecture = computer('arch');
metadata.matlab.toolboxes = ver;

metadata.input.file = 'tests/test_data.csv';
metadata.input.sha256 = file_sha256(input_path);
metadata.input.row_count = height(input_data);
metadata.input.valid_isotope_pair_count = sum(valid_rows);
metadata.input.variables = input_data.Properties.VariableNames;

metadata.configuration = config;
metadata.interpolation.temperature_method = 'natural';
metadata.interpolation.mixing_ratio_method = 'linear (MATLAB griddata default)';
metadata.randomness = 'none';

active_functions = {
    'simple_water_isotope_model_2020'
    'evaporation_2021'
    'T_RH_RHn_2020'
    'd18Osw_to_dDsw'
    'distillation_2020'
    'fraction_il_brm_H10'
    'mixed_phased_supersaturation'
    'pseudo_adiabat_function'
    'Ts_to_Tc_2020'
    };

source_files = repmat(struct('function_name', '', 'file', '', 'sha256', ''), ...
    numel(active_functions) + 1, 1);
source_files(1).function_name = 'run_swim_allan_hills_test';
source_files(1).file = 'tests/run_swim_allan_hills_test.m';
source_files(1).sha256 = file_sha256(script_path);

for source_index = 1:numel(active_functions)
    function_name = active_functions{source_index};
    function_path = which(function_name);
    assert(~isempty(function_path), 'Required function is not on the path: %s', ...
        function_name);
    source_files(source_index + 1).function_name = function_name;
    source_files(source_index + 1).file = ...
        ['legacy_matlab/' function_name '.m'];
    source_files(source_index + 1).sha256 = file_sha256(function_path);
end
metadata.source_files = source_files;

%% Run the legacy forward model

% T_RH_RHn_2020 loads climatology files using paths relative to
% legacy_matlab, so run the model from that directory and restore the
% caller's working directory afterward.
original_directory = pwd;
directory_cleanup = onCleanup(@() cd(original_directory));
cd(legacy_dir);

run_timer = tic;
[T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, ...
    dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site, P_site] = ...
    simple_water_isotope_model_2020( ...
        config.Tsite, ...
        config.Tsource, ...
        config.RHsource, ...
        config.supersaturation_a, ...
        config.supersaturation_b, ...
        config.supersaturation_c, ...
        config.closure, ...
        config.reanalysis, ...
        config.hemisphere, ...
        config.season);
model_runtime_seconds = toc(run_timer);

%% Run one detailed diagnostic trajectory

trajectory_timer = tic;
[trajectory_dD_v0, trajectory_d18O_v0, trajectory_d17Oxs_v0, ...
    trajectory_RHn0, trajectory_RH0, trajectory_SST0] = ...
    evaporation_2021( ...
        config.trajectory.Tsource_C, [], [], config.closure, ...
        config.reanalysis, config.hemisphere, config.season);

[trajectory_d18O_p, trajectory_dD_p, trajectory_d17O_p, ...
    trajectory_d18O_p_ln, trajectory_dD_p_ln, trajectory_d17O_p_ln, ...
    trajectory_dxs, trajectory_dln, trajectory_d17Oxs, trajectory_T, ...
    trajectory_S, trajectory_supersaturation, trajectory_r_s, ...
    trajectory_e_s, trajectory_f, trajectory_P, trajectory_RH, ...
    trajectory_RD_vapor, trajectory_RD_precip, ...
    trajectory_R18O_vapor, trajectory_R18O_precip, ...
    trajectory_R17O_vapor, trajectory_R17O_precip] = ...
    distillation_2020( ...
        config.trajectory.Tsource_C, config.trajectory.Tcond_C, ...
        config.trajectory.step_degC, trajectory_dD_v0, trajectory_d18O_v0, ...
        trajectory_d17Oxs_v0, trajectory_RHn0, trajectory_RH0, ...
        trajectory_SST0, config.supersaturation_a, ...
        config.supersaturation_b, config.supersaturation_c, 'adiabatic');

[trajectory_fraction_ice, trajectory_fraction_liquid] = ...
    fraction_il_brm_H10(trajectory_T, 'adj');
trajectory_alpha_D_effective = trajectory_RD_precip ./ trajectory_RD_vapor;
trajectory_alpha_18O_effective = ...
    trajectory_R18O_precip ./ trajectory_R18O_vapor;
trajectory_alpha_17O_effective = ...
    trajectory_R17O_precip ./ trajectory_R17O_vapor;
trajectory_runtime_seconds = toc(trajectory_timer);

cd(original_directory);
clear directory_cleanup;

%% Reconstruct temperatures from the Allan Hills observations

observed_d18O_ln = log(1 + observed_d18O ./ 1000) .* 1000;
observed_dD_ln = log(1 + observed_dD ./ 1000) .* 1000;
observed_dln = observed_dD_ln - ...
    (-2.85e-2 .* observed_d18O_ln .^ 2 + 8.47 .* observed_d18O_ln);
observed_dxs_calculated = observed_dD - 8 .* observed_d18O;

source_grid = repmat(T_source, length(T_site), 1)';
site_grid = repmat(T_site, length(T_source), 1);
state_valid = ~isnan(d18Oln_site);

reconstructed_Tcond_C = nan(height(input_data), 1);
reconstructed_Tsource_C = nan(height(input_data), 1);
reconstructed_Tsurface_C = nan(height(input_data), 1);
reconstructed_r_s = nan(height(input_data), 1);

reconstructed_Tcond_C(valid_rows) = griddata( ...
    d18Oln_site(state_valid), dlnU_site(state_valid), site_grid(state_valid), ...
    observed_d18O_ln(valid_rows), observed_dln(valid_rows), 'natural');
reconstructed_Tsource_C(valid_rows) = griddata( ...
    d18Oln_site(state_valid), dlnU_site(state_valid), source_grid(state_valid), ...
    observed_d18O_ln(valid_rows), observed_dln(valid_rows), 'natural');
reconstructed_r_s(valid_rows) = griddata( ...
    d18Oln_site(state_valid), dlnU_site(state_valid), r_s_site(state_valid), ...
    observed_d18O_ln(valid_rows), observed_dln(valid_rows));
reconstructed_Tsurface_C(valid_rows) = ...
    Ts_to_Tc_2020(reconstructed_Tcond_C(valid_rows), 1);

results = input_data;
results.d18O_ln = observed_d18O_ln;
results.dD_ln = observed_dD_ln;
results.dln = observed_dln;
results.dxs_calculated = observed_dxs_calculated;
results.Tcond_reconstructed_C = reconstructed_Tcond_C;
results.Tsource_reconstructed_C = reconstructed_Tsource_C;
results.Tsurface_reconstructed_C = reconstructed_Tsurface_C;
results.r_s_reconstructed = reconstructed_r_s;

if ismember('dxs_smow', input_data.Properties.VariableNames)
    results.dxs_residual = input_data.dxs_smow - observed_dxs_calculated;
    metadata.input.max_abs_dxs_residual = max(abs(results.dxs_residual), ...
        [], 'omitnan');
end

metadata.result.reconstructed_count = sum( ...
    isfinite(reconstructed_Tcond_C) & isfinite(reconstructed_Tsource_C));
metadata.result.outside_interpolation_domain_count = sum( ...
    valid_rows & (~isfinite(reconstructed_Tcond_C) | ...
    ~isfinite(reconstructed_Tsource_C)));
metadata.result.invalid_input_row_count = sum(~valid_rows);
metadata.result.model_runtime_seconds = model_runtime_seconds;
metadata.result.state_space_shape = size(d18O_site);
metadata.result.trajectory_runtime_seconds = trajectory_runtime_seconds;
metadata.result.trajectory_row_count = numel(trajectory_T);
metadata.result.trajectory_initial_conditions.dD_v0 = trajectory_dD_v0;
metadata.result.trajectory_initial_conditions.d18O_v0 = trajectory_d18O_v0;
metadata.result.trajectory_initial_conditions.d17Oxs_v0 = trajectory_d17Oxs_v0;
metadata.result.trajectory_initial_conditions.RHn0 = trajectory_RHn0;
metadata.result.trajectory_initial_conditions.RH0 = trajectory_RH0;
metadata.result.trajectory_initial_conditions.SST0_C = trajectory_SST0;

finished_at = datetime('now', 'TimeZone', 'UTC');
finished_at.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX';
metadata.finished_at_utc = char(finished_at);
metadata.status = 'completed';

%% Save outputs

results_csv_path = fullfile(output_dir, 'allan_hills_reconstruction.csv');
state_space_mat_path = fullfile(output_dir, 'state_space.mat');
trajectory_csv_path = fullfile(output_dir, 'trajectory_Tsource_10_Tcond_m30.csv');
trajectory_mat_path = fullfile(output_dir, 'trajectory_Tsource_10_Tcond_m30.mat');
metadata_json_path = fullfile(output_dir, 'metadata.json');

writetable(results, results_csv_path);

% Keep MATLAB tables and structs out of these .mat files so scipy.io.loadmat
% can inspect and load them without MATLAB-specific object decoding.
save(state_space_mat_path, ...
    'T_site', 'T_source', 'RH_source', ...
    'd18O_site', 'dD_site', 'd18Oln_site', 'dDln_site', ...
    'dxs_site', 'd17O_xs_site', 'dlnU_site', 'r_s_site', 'P_site', ...
    '-v7');

trajectory = table;
trajectory.temperature_C = trajectory_T(:);
trajectory.pressure_kPa = trajectory_P(:);
trajectory.saturation_vapor_pressure_kPa = trajectory_e_s(:);
trajectory.saturated_mixing_ratio = trajectory_r_s(:);
trajectory.fraction_vapor_remaining = trajectory_f(:);
trajectory.saturation_used = trajectory_S(:);
trajectory.mixed_phase_supersaturation = trajectory_supersaturation(:);
trajectory.relative_humidity = trajectory_RH(:);
trajectory.fraction_ice = trajectory_fraction_ice(:);
trajectory.fraction_liquid = trajectory_fraction_liquid(:);
trajectory.d18O_precip = trajectory_d18O_p(:);
trajectory.dD_precip = trajectory_dD_p(:);
trajectory.d17O_precip = trajectory_d17O_p(:);
trajectory.d18O_ln = trajectory_d18O_p_ln(:);
trajectory.dD_ln = trajectory_dD_p_ln(:);
trajectory.d17O_ln = trajectory_d17O_p_ln(:);
trajectory.dxs = trajectory_dxs(:);
trajectory.dln = trajectory_dln(:);
trajectory.d17O_excess = trajectory_d17Oxs(:);
trajectory.RD_vapor = trajectory_RD_vapor(:);
trajectory.RD_precip = trajectory_RD_precip(:);
trajectory.R18O_vapor = trajectory_R18O_vapor(:);
trajectory.R18O_precip = trajectory_R18O_precip(:);
trajectory.R17O_vapor = trajectory_R17O_vapor(:);
trajectory.R17O_precip = trajectory_R17O_precip(:);
trajectory.alpha_D_effective = trajectory_alpha_D_effective(:);
trajectory.alpha_18O_effective = trajectory_alpha_18O_effective(:);
trajectory.alpha_17O_effective = trajectory_alpha_17O_effective(:);
writetable(trajectory, trajectory_csv_path);

save(trajectory_mat_path, ...
    'trajectory_T', 'trajectory_P', 'trajectory_e_s', 'trajectory_r_s', ...
    'trajectory_f', 'trajectory_S', 'trajectory_supersaturation', ...
    'trajectory_RH', 'trajectory_fraction_ice', ...
    'trajectory_fraction_liquid', 'trajectory_d18O_p', ...
    'trajectory_dD_p', 'trajectory_d17O_p', 'trajectory_d18O_p_ln', ...
    'trajectory_dD_p_ln', 'trajectory_d17O_p_ln', 'trajectory_dxs', ...
    'trajectory_dln', 'trajectory_d17Oxs', 'trajectory_RD_vapor', ...
    'trajectory_RD_precip', 'trajectory_R18O_vapor', ...
    'trajectory_R18O_precip', 'trajectory_R17O_vapor', ...
    'trajectory_R17O_precip', 'trajectory_alpha_D_effective', ...
    'trajectory_alpha_18O_effective', 'trajectory_alpha_17O_effective', ...
    'trajectory_dD_v0', 'trajectory_d18O_v0', 'trajectory_d17Oxs_v0', ...
    'trajectory_RHn0', 'trajectory_RH0', 'trajectory_SST0', '-v7');

metadata.outputs.results_csv.file = 'allan_hills_reconstruction.csv';
metadata.outputs.results_csv.sha256 = file_sha256(results_csv_path);
metadata.outputs.state_space_mat.file = 'state_space.mat';
metadata.outputs.state_space_mat.sha256 = file_sha256(state_space_mat_path);
metadata.outputs.trajectory_csv.file = 'trajectory_Tsource_10_Tcond_m30.csv';
metadata.outputs.trajectory_csv.sha256 = file_sha256(trajectory_csv_path);
metadata.outputs.trajectory_mat.file = 'trajectory_Tsource_10_Tcond_m30.mat';
metadata.outputs.trajectory_mat.sha256 = file_sha256(trajectory_mat_path);

metadata_json = jsonencode(metadata, 'PrettyPrint', true);
write_text_file(metadata_json_path, metadata_json);

fprintf('Completed in %.1f seconds.\n', model_runtime_seconds);
fprintf('Reconstructed %d of %d valid isotope pairs.\n', ...
    metadata.result.reconstructed_count, metadata.input.valid_isotope_pair_count);
fprintf('Results: %s\n', results_csv_path);
fprintf('State space: %s\n', state_space_mat_path);
fprintf('Trajectory CSV: %s\n', trajectory_csv_path);
fprintf('Trajectory MAT: %s\n', trajectory_mat_path);
fprintf('Metadata: %s\n', metadata_json_path);

%% Local helpers

function value = git_value(project_root, arguments)
    escaped_root = strrep(project_root, '"', '\"');
    command = sprintf('git -C "%s" %s', escaped_root, arguments);
    [status, output] = system(command);
    if status == 0
        value = strtrim(output);
    else
        value = '';
    end
end

function checksum = file_sha256(file_path)
    file_id = fopen(file_path, 'rb');
    assert(file_id >= 0, 'Could not open file for hashing: %s', file_path);
    file_cleanup = onCleanup(@() fclose(file_id));
    bytes = fread(file_id, Inf, '*uint8');

    digest = java.security.MessageDigest.getInstance('SHA-256');
    digest.update(bytes);
    checksum_bytes = typecast(digest.digest(), 'uint8');
    checksum = lower(reshape(dec2hex(checksum_bytes, 2)', 1, []));

    clear file_cleanup;
end

function write_text_file(file_path, contents)
    file_id = fopen(file_path, 'w');
    assert(file_id >= 0, 'Could not open text file for writing: %s', file_path);
    file_cleanup = onCleanup(@() fclose(file_id));
    fprintf(file_id, '%s\n', contents);
    clear file_cleanup;
end
