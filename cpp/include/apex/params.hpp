// GENERE par scripts/gen_params.py depuis les dataclasses Python. NE PAS EDITER.
// Contrat de parametres du moteur (meme contenu que config/params.schema.json).
// Les champs Python `typing.Any` (ex. MachineProfile.nozzle_model) sont
// volontairement absents : donnees d'execution, pas parametres.
#pragma once

#include <array>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace apex::params {

enum class ConstraintFalloff { LINEAR, GAUSSIAN, CUBIC };
inline const char* to_string(ConstraintFalloff v) {
    switch (v) {
        case ConstraintFalloff::LINEAR: return "linear";
        case ConstraintFalloff::GAUSSIAN: return "gaussian";
        case ConstraintFalloff::CUBIC: return "cubic";
    }
    return "";
}

enum class ConstraintType { ATTACH_POINT, AESTHETIC, STRUCTURAL };
inline const char* to_string(ConstraintType v) {
    switch (v) {
        case ConstraintType::ATTACH_POINT: return "attach_point";
        case ConstraintType::AESTHETIC: return "aesthetic";
        case ConstraintType::STRUCTURAL: return "structural";
    }
    return "";
}

enum class InfillPattern { LINES, GRID, CROSSED_ZIGZAG, OSCILLATING_SPIRAL };
inline const char* to_string(InfillPattern v) {
    switch (v) {
        case InfillPattern::LINES: return "lines";
        case InfillPattern::GRID: return "grid";
        case InfillPattern::CROSSED_ZIGZAG: return "crossed_zigzag";
        case InfillPattern::OSCILLATING_SPIRAL: return "oscillating_spiral";
    }
    return "";
}

struct OverhangParams {
    double overhang_angle_threshold_deg{45.0};
    int bridge_detection_ray_density{4};
    double bridge_support_ratio_threshold{0.7};
    double bridge_max_gap_distance{2.0};
    double bridge_horizontal_angle_tolerance_deg{10.0};
    std::vector<double> up_vector{std::vector<double>{0.0, 0.0, 1.0}};
};

struct RotationAxisSpec {
    std::vector<double> axis{};
    std::array<double, 2> angle_range_deg{};
};

struct SegmentationParams {
    int direction_candidate_count{26};
    double max_angle_from_direction_deg{45.0};
    std::optional<double> max_region_merge_angle_deg{std::nullopt};
    std::optional<double> platform_z{std::nullopt};
    double platform_contact_tolerance{0.5};
    std::optional<std::vector<RotationAxisSpec>> rotation_axes{std::nullopt};
    std::vector<double> up_vector{std::vector<double>{0.0, 0.0, 1.0}};
    std::optional<double> anchor_growth_weight{std::nullopt};
    std::optional<double> branch_cut_threshold{std::nullopt};
};

struct FieldSolverParams {
    double laplacian_smoothness_weight{1.0};
    double gradient_overhang_weight{5.0};
    double resolution{1.5};
    int max_solver_iterations{2000};
    double convergence_tolerance{1e-08};
    double gradient_baseline_weight{1.0};
    int overhang_influence_depth{1};
    bool enable_platform_anchor{true};
    double platform_anchor_weight{1000000.0};
    double platform_anchor_value{0.0};
    std::optional<std::vector<double>> up_vector{std::vector<double>{0.0, 0.0, 1.0}};
    int direct_solve_node_threshold{30000};
    bool use_ilu_preconditioner{false};
    double direction_blend_radius{0.0};
    double direction_blend_power{2.0};
    double direction_blend_min_norm{0.1};
    bool enforce_upward_growth{false};
    double upward_growth_epsilon{1e-09};
    std::optional<double> max_downward_growth_angle_deg{std::nullopt};
    bool use_geodesic_directions{false};
    double geodesic_source_weight{0.0};
    double geodesic_saddle_threshold_deg{100.0};
    int geodesic_smoothing_iterations{0};
    double geodesic_smoothing_step{0.5};
    int support_feedback_iterations{0};
    std::optional<std::vector<double>> radial_tilt_origin{std::nullopt};
    double support_feedback_min_improvement{0.02};
    double support_feedback_threshold_deg{45.0};
    double support_feedback_weight{4.0};
    bool support_feedback_use_geodesic_target{false};
    double smoothness_modulation_strength{0.0};
    double smoothness_modulation_floor{0.02};
    std::optional<double> transition_distance_mm{std::nullopt};
    double heat_method_time_factor{100.0};
    bool use_heat_geodesic_field{false};
    bool heat_project_reachable{false};
    bool multibranch{false};
    bool planar_parts{false};
    std::string planar_parts_partition{"region"};
    double planar_reeb_band_mm{4.0};
    int multibranch_min_branch_tets{40};
    std::string multibranch_direction_mode{"region"};
    int multibranch_smoothing_iterations{0};
    std::string multibranch_partition{"region"};
    double multibranch_flow_angle_deg{30.0};
    bool multibranch_use_risk_weights{true};
    bool multibranch_order_by_branch{true};
    bool multibranch_blend_directions{false};
    std::optional<double> multibranch_vertical_if_safe{std::nullopt};
    bool support_feedback_require_in_air{false};
    double support_feedback_in_air_threshold_mm{1.0};
};

struct Constraint {
    ConstraintType type{};
    std::vector<double> position{};
    std::optional<double> target_field_value{};
    std::optional<std::vector<double>> target_normal{};
    double influence_radius{};
    ConstraintFalloff falloff{};
    double weight{};
};

struct LayerParams {
    double layer_height_base{0.2};
    double weld_tolerance_factor{1e-06};
    std::optional<double> split_disjoint_clusters{std::nullopt};
    std::optional<double> support_driven_order{std::nullopt};
    bool geodesic_order{false};
    bool frontier_order{false};
    bool isovalue_order{false};
    std::optional<double> nozzle_width{std::nullopt};
    double support_horizontal_radius_factor{2.5};
    double support_vertical_height_factor{7.5};
    double platform_contact_tolerance_factor{2.5};
    std::optional<std::vector<RotationAxisSpec>> support_rotation_axes{std::nullopt};
    int support_direction_candidate_count{26};
};

struct Step6Params {
    double nozzle_width{0.4};
    int wall_count{2};
    double infill_density{0.2};
    InfillPattern infill_pattern{InfillPattern::CROSSED_ZIGZAG};
    double bridge_contour_ratio_threshold{0.05};
    int bridge_anchor_extra_walls{2};
    double degenerate_area_ratio{1e-06};
    std::optional<double> out_of_plane_interpolation_power{std::nullopt};
    bool split_rings_into_runs{false};
    bool split_open_polylines{false};
    bool continuous_rings{false};
    double run_split_min_amplitude{1e-12};
    std::optional<std::array<double, 3>> spiral_axis_origin{std::nullopt};
    std::optional<double> spiral_alternation_threshold{std::nullopt};
    std::optional<double> spiral_peak_spacing{std::nullopt};
    double spiral_randomization_threshold{0.5};
    int spiral_seed{0};
};

struct HoleCostModel {
    double reheat_penalty{1.0};
    double adhesion_risk_weight{1.0};
};

struct PurgeCostModel {
    double volume_per_purge{1.0};
    double time_per_purge{1.0};
};

struct Step7Params {
    double surface_color_distance{0.5};
    HoleCostModel hole_cost_model{HoleCostModel{.reheat_penalty = 1.0, .adhesion_risk_weight = 1.0}};
    PurgeCostModel purge_cost_model{PurgeCostModel{.volume_per_purge = 1.0, .time_per_purge = 1.0}};
    int max_open_holes{4};
    int escape_direction_count{16};
    std::optional<std::vector<RotationAxisSpec>> rotation_axes{std::nullopt};
    int deadline_ray_samples{24};
};

struct ColorMap {
    std::map<int, int> per_face_color{};
    double merge_tolerance{};
};

struct Step8Params {
    int ik_samples_per_axis{10};
    int ik_refine_iterations{8};
    double ik_refine_step_multiplier{0.5};
    double singularity_radius_epsilon{0.001};
    double travel_feedrate{0.0};
    double default_feedrate{1000.0};
    double max_unreachable_angular_error_deg{5.0};
    bool clamp_angle_rates{false};
    int free_axis_median_window{9};
    bool per_point_normal{true};
    bool per_point_normal_nodal{false};
};

struct CollisionSimParams {
    double path_sampling_step{};
    double safety_margin{};
    int max_adjustment_iterations{};
    double grid_cell_size_multiplier{};
    double toolhead_sample_spacing{};
    double nozzle_contact_clearance{};
    bool fast_kdtree{};
    double fast_densify_step{};
};

struct Step10Params {
    int max_adjustment_iterations{3};
    int reorientation_samples_per_axis{12};
    double max_direction_deviation_deg{45.0};
    double excluded_range_half_width_deg{15.0};
    int rate_limit_neighbor_margin{1};
    std::optional<double> orientation_support_min_residual_deg{std::nullopt};
};

struct AABB {
    std::vector<double> min_corner{};
    std::vector<double> max_corner{};
};

struct MachineProfile {
    std::vector<RotationAxisSpec> rotation_axes{std::vector<RotationAxisSpec>{}};
    std::vector<double> platform_origin{};
    // nozzle_model : donnee d'execution, non transcrite
    double nozzle_diameter{0.4};
    double max_angular_velocity{0.0};
    double max_angular_acceleration{0.0};
    std::optional<AABB> build_volume{std::nullopt};
    double z_clearance_min{0.0};
    std::optional<double> platform_radius{std::nullopt};
    std::vector<int> position_coupled_axes{std::vector<int>{}};
    std::optional<std::vector<double>> radial_direction{std::nullopt};
};

struct TravelRoutingParams {
    double length_threshold{5.0};
    double lift_height{3.0};
    double max_lift_multiplier{4.0};
    int traverse_samples{1};
    double routing_toolhead_sample_spacing{15.0};
    double routing_grid_cell_size_multiplier{40.0};
};

struct CollisionOrderParams {
    double pose_spacing{3.0};
    double material_spacing{1.5};
    double head_sample_spacing{6.0};
    double margin{1.0};
    double nozzle_contact_clearance{3.0};
    int neighbours{4};
    std::optional<double> support_depth{0.3};
    double support_radius{0.6};
};

struct SupportCheckParams {
    double horizontal_search_radius_factor{2.5};
    double vertical_search_height_factor{7.5};
    double platform_contact_tolerance_factor{2.5};
    double path_sampling_step{0.0};
    bool locate_supporting_point{false};
    double grid_cell_size_multiplier{1.0};
    int direction_candidate_count{26};
};

struct PipelineParams {
    double resolution{1.5};
    bool allow_multi_mesh{false};
    std::optional<double> near_platform_resolution{0.5};
    double near_platform_margin{5.0};
    std::optional<double> junction_refinement_resolution{std::nullopt};
    double junction_refinement_angle_deg{30.0};
    double junction_refinement_margin{3.0};
    bool place_on_platform{true};
    OverhangParams overhang{OverhangParams{.overhang_angle_threshold_deg = 45.0, .bridge_detection_ray_density = 4, .bridge_support_ratio_threshold = 0.7, .bridge_max_gap_distance = 2.0, .bridge_horizontal_angle_tolerance_deg = 10.0, .up_vector = std::vector<double>{0.0, 0.0, 1.0}}};
    SegmentationParams segmentation{SegmentationParams{.direction_candidate_count = 26, .max_angle_from_direction_deg = 45.0, .max_region_merge_angle_deg = std::nullopt, .platform_z = std::nullopt, .platform_contact_tolerance = 0.5, .rotation_axes = std::nullopt, .up_vector = std::vector<double>{0.0, 0.0, 1.0}, .anchor_growth_weight = std::nullopt, .branch_cut_threshold = std::nullopt}};
    FieldSolverParams field_solver{FieldSolverParams{.laplacian_smoothness_weight = 1.0, .gradient_overhang_weight = 5.0, .resolution = 1.5, .max_solver_iterations = 2000, .convergence_tolerance = 1e-08, .gradient_baseline_weight = 1.0, .overhang_influence_depth = 1, .enable_platform_anchor = true, .platform_anchor_weight = 1000000.0, .platform_anchor_value = 0.0, .up_vector = std::vector<double>{0.0, 0.0, 1.0}, .direct_solve_node_threshold = 30000, .use_ilu_preconditioner = false, .direction_blend_radius = 10.0, .direction_blend_power = 2.0, .direction_blend_min_norm = 0.1, .enforce_upward_growth = true, .upward_growth_epsilon = 1e-09, .max_downward_growth_angle_deg = std::nullopt, .use_geodesic_directions = false, .geodesic_source_weight = 0.0, .geodesic_saddle_threshold_deg = 100.0, .geodesic_smoothing_iterations = 0, .geodesic_smoothing_step = 0.5, .support_feedback_iterations = 0, .radial_tilt_origin = std::nullopt, .support_feedback_min_improvement = 0.02, .support_feedback_threshold_deg = 45.0, .support_feedback_weight = 4.0, .support_feedback_use_geodesic_target = false, .smoothness_modulation_strength = 0.0, .smoothness_modulation_floor = 0.02, .transition_distance_mm = std::nullopt, .heat_method_time_factor = 100.0, .use_heat_geodesic_field = false, .heat_project_reachable = false, .multibranch = false, .planar_parts = false, .planar_parts_partition = "region", .planar_reeb_band_mm = 4.0, .multibranch_min_branch_tets = 40, .multibranch_direction_mode = "region", .multibranch_smoothing_iterations = 0, .multibranch_partition = "region", .multibranch_flow_angle_deg = 30.0, .multibranch_use_risk_weights = true, .multibranch_order_by_branch = true, .multibranch_blend_directions = false, .multibranch_vertical_if_safe = std::nullopt, .support_feedback_require_in_air = false, .support_feedback_in_air_threshold_mm = 1.0}};
    std::optional<std::vector<Constraint>> constraints{std::nullopt};
    LayerParams layers{LayerParams{.layer_height_base = 0.2, .weld_tolerance_factor = 1e-06, .split_disjoint_clusters = std::nullopt, .support_driven_order = std::nullopt, .geodesic_order = false, .frontier_order = false, .isovalue_order = false, .nozzle_width = std::nullopt, .support_horizontal_radius_factor = 2.5, .support_vertical_height_factor = 7.5, .platform_contact_tolerance_factor = 2.5, .support_rotation_axes = std::nullopt, .support_direction_candidate_count = 26}};
    Step6Params toolpath2d{Step6Params{.nozzle_width = 0.4, .wall_count = 2, .infill_density = 0.2, .infill_pattern = InfillPattern::CROSSED_ZIGZAG, .bridge_contour_ratio_threshold = 0.05, .bridge_anchor_extra_walls = 2, .degenerate_area_ratio = 1e-06, .out_of_plane_interpolation_power = 2.0, .split_rings_into_runs = true, .split_open_polylines = false, .continuous_rings = true, .run_split_min_amplitude = 0.05, .spiral_axis_origin = std::nullopt, .spiral_alternation_threshold = std::nullopt, .spiral_peak_spacing = std::nullopt, .spiral_randomization_threshold = 0.5, .spiral_seed = 0}};
    Step7Params scheduling{Step7Params{.surface_color_distance = 0.5, .hole_cost_model = HoleCostModel{.reheat_penalty = 1.0, .adhesion_risk_weight = 1.0}, .purge_cost_model = PurgeCostModel{.volume_per_purge = 1.0, .time_per_purge = 1.0}, .max_open_holes = 4, .escape_direction_count = 16, .rotation_axes = std::nullopt, .deadline_ray_samples = 24}};
    std::optional<ColorMap> color_map{std::nullopt};
    Step8Params machine_toolpath{Step8Params{.ik_samples_per_axis = 10, .ik_refine_iterations = 8, .ik_refine_step_multiplier = 0.5, .singularity_radius_epsilon = 0.001, .travel_feedrate = 0.0, .default_feedrate = 1000.0, .max_unreachable_angular_error_deg = 5.0, .clamp_angle_rates = false, .free_axis_median_window = 9, .per_point_normal = true, .per_point_normal_nodal = false}};
    CollisionSimParams collision{CollisionSimParams{.path_sampling_step = 2.0, .safety_margin = 1.0, .max_adjustment_iterations = 3, .grid_cell_size_multiplier = 4.0, .toolhead_sample_spacing = 5.0, .nozzle_contact_clearance = 3.0, .fast_kdtree = true, .fast_densify_step = 0.25}};
    Step10Params adjustment{Step10Params{.max_adjustment_iterations = 3, .reorientation_samples_per_axis = 12, .max_direction_deviation_deg = 45.0, .excluded_range_half_width_deg = 15.0, .rate_limit_neighbor_margin = 1, .orientation_support_min_residual_deg = std::nullopt}};
    bool height_clamp{false};
    std::optional<double> firmware_z_min{std::nullopt};
    std::optional<std::vector<double>> coupled_reorientation_tiers{std::nullopt};
    int coupled_reorientation_rounds{3};
    MachineProfile machine_profile{MachineProfile{.rotation_axes = std::vector<RotationAxisSpec>{RotationAxisSpec{.axis = std::vector<double>{0.0, 0.0, 1.0}, .angle_range_deg = std::array<double, 2>{0.0, 360.0}}, RotationAxisSpec{.axis = std::vector<double>{1.0, 0.0, 0.0}, .angle_range_deg = std::array<double, 2>{0.0, 185.0}}}, .platform_origin = std::vector<double>{0.0, 0.0, 0.0}, .nozzle_diameter = 0.4, .max_angular_velocity = 29.08940264298949, .max_angular_acceleration = 22.689280275926283, .build_volume = std::nullopt, .z_clearance_min = 0.0, .platform_radius = 150.0, .position_coupled_axes = std::vector<int>{0}, .radial_direction = std::vector<double>{0.0, -1.0, 0.0}}};
    std::optional<TravelRoutingParams> travel_routing{std::nullopt};
    std::optional<CollisionOrderParams> collision_order{CollisionOrderParams{.pose_spacing = 3.0, .material_spacing = 1.5, .head_sample_spacing = 6.0, .margin = 1.0, .nozzle_contact_clearance = 3.0, .neighbours = 4, .support_depth = 0.3, .support_radius = 0.6}};
    std::optional<SupportCheckParams> support_check{std::nullopt};
    bool reachable_tilt_from_machine{true};
};

struct GcodeParams {
    double nozzle_temp{195.0};
    double filament_diameter{1.75};
    double layer_height{0.2};
    double skin_thickness{0.2};
    double line_width{0.4};
    double flow{1.0};
    double print_speed{25.0};
    double travel_speed{80.0};
    double rotation_feed{3000.0};
    double max_feed{12000.0};
    double retract_length{0.8};
    double retract_speed{35.0};
    double retract_min_travel{2.0};
    double travel_clearance{1.0};
    double travel_clearance_max_factor{4.0};
    double pivot_to_tip{46.0};
    double tilt_sign{1.0};
    std::string z_reference{"tip"};
    double b_sign{1.0};
    bool tip_at_bead_top{true};
    bool strict{false};
    std::optional<double> max_y_turns{std::nullopt};
    double platform_margin{1.0};
    double singularity_radius{1e-06};
    int decimals{4};
    std::optional<std::vector<double>> fw_axis_speeds{std::vector<double>{1666.7, 600.0, 133.3, 1666.7}};
    std::optional<double> max_tip_deviation{0.05};
};

}  // namespace apex::params
