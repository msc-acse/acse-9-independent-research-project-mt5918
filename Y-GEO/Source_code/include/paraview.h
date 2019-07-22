/*

  Paraview output for FEMDEM.

*/


#ifndef paraview_H
#define paraview_H

/* Paraview vector */
typedef struct 
{
    float x;
    float y;
    float z;
} pv_vector_t;

/* Paraview tensor */
typedef struct
{
    float a11, a12, a13, a21, a22, a23, a31, a32, a33;
} pv_tensor_t;


/* Initialization function - Must be called before any other calls to pv_ functions */
void pv_init(int num_elements, int upper_bound_joints, int num_rebars, int upper_bound_1D_joints, int num_dead_nodes, int upper_bound_acoustic_events, const char* file_prefix, double realtime_per_timestep);

/* Destruction function - Must be called after all calls to pv_ functions are done */
void pv_free();

/* Write output file for the timestep.  Called when the timestep is done. */
void pv_write_timestep(int timestep);

/* Data output functions */

/* element_id and joint_id are munjiza mangled */
void pv_write_stress(int munjiza_element_id, pv_tensor_t tensor); 
void pv_write_principal_stress(int munjiza_element_id, pv_tensor_t tensor); 
void pv_write_strain(int munjiza_element_id, pv_tensor_t tensor); 
//void pv_write_principal_stress_direction_1(int munjiza_element_id, pv_vector_t vector);
//void pv_write_principal_stress_direction_2(int munjiza_element_id, pv_vector_t vector);

void pv_write_nodal_force(int node_id, pv_vector_t vector);
void pv_write_nodal_current_coord(int node_id, pv_vector_t vector);
void pv_write_nodal_velocity(int node_id, pv_vector_t vector);
void pv_write_nodal_displacement(int node_id, pv_vector_t vector);
void pv_write_nodal_fluid_pressure(int node_id, float fluid_pressure);

void pv_write_element_property_id(int munjiza_element_id, int property_id);

void pv_write_broken_joint_mode_of_failure(int broken_joint_id, float mode_of_failure);
void pv_write_broken_joint_sliding(int broken_joint_id, float mode_of_failure);
void pv_write_broken_joint_opening(int broken_joint_id, float mode_of_failure);
void pv_write_broken_joint_area(int broken_joint_id, float mode_of_failure);
void pv_write_yielded_joint_mode_of_failure(int yielded_joint_id, int mode_of_failure);

void pv_write_broken_joint_coordinates(int broken_joint_id, pv_vector_t pt0,  pv_vector_t pt1,  pv_vector_t pt2,  pv_vector_t pt3); 
void pv_write_yielded_joint_coordinates(int yielded_joint_id, pv_vector_t pt0,  pv_vector_t pt1,  pv_vector_t pt2,  pv_vector_t pt3); 

void pv_write_rebar_coordinates(int broken_joint_id, pv_vector_t pt0,  pv_vector_t pt1); 
void pv_write_rebar_area(int rebar_id, float area);
void pv_write_rebar_stress(int rebar_id, float stress);
void pv_write_rebar_force(int rebar_id, float force);
void pv_write_rebar_strain(int rebar_id, float strain);
void pv_write_rebar_activation_flag(int rebar_id, int rebar_activation_flag);


void pv_write_1D_joint_coordinates(int oneD_joint_id, pv_vector_t pt0,  pv_vector_t pt1);
void pv_write_1D_joint_normal_stress(int oneD_joint_id, float stress);
void pv_write_1D_joint_tangential_stress(int oneD_joint_id, float stress);
void pv_write_1D_joint_normal_strain(int oneD_joint_id, float strain);
void pv_write_1D_joint_tangential_strain(int oneD_joint_id, float strain);
void pv_write_1D_joint_state(int oneD_joint_id, int state);
//void pv_write_1D_joint_force(int oneD_joint_id, pv_vector_t vector);

void pv_set_num_broken_joints(int count);
void pv_set_num_yielded_joints(int count);
void pv_set_num_acoustic_events(int count);
void pv_set_num_unbroken_1D_joints(int count);

void pv_write_element_connectivity(int munjiza_element_id, int node1, int node2, int node3);

/* Acoustic emissions */
void pv_write_acoustic_coordinate(int acoustic_node_id, pv_vector_t coord);
void pv_write_acoustic_energy(int acoustic_node_id, float energy);
void pv_write_acoustic_magnitude(int acoustic_node_id, float magnitude);
void pv_write_acoustic_fracture_energy(int acoustic_node_id, float energy);
void pv_write_acoustic_event_mode(int acoustic_node_id, float mode);

/* Principal directions */
void pv_write_element_centroid_coordinate(int munjiza_element_id, pv_vector_t coord);
void pv_write_principal_stress_direction_1_andrea(int munjiza_element_id, pv_vector_t vector);
void pv_write_principal_stress_direction_2_andrea(int munjiza_element_id, pv_vector_t vector);


/* De-Munjiza Mangler */

/* For each joint in the simulation, let us know how munjiza labels it */
void pv_register_munjiza_joint(int munjiza_joint);

/* For each element in the simulation, let us know how munjiza labels it */
void pv_register_munjiza_element(int munjiza_element);

/* Translate Munjiza element ID to sane ID */
int translate_munjiza_element_to_sane(int munjiza_element);

/* Translate Munjiza joint ID to sane ID */
int translate_munjiza_joint_to_sane(int munjiza_joint);

#endif
