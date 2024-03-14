/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************/

#ifndef FEA_MODULE_SGH_H
#define FEA_MODULE_SGH_H

#include "mesh.h"
#include "state.h"
#include "matar.h"
// #include "elements.h"
#include "solver.h"
// #include "FEA_Module.h"
// #include "material_models.h"



struct material_t;

struct boundary_t;

/////////////////////////////////////////////////////////////////////////////
///
/// \class SGH
///
/// \brief Class for containing functions required to perform SGH
///
/// This class containts the requisite functions requited to perform
/// staggered grid hydrodynamics (SGH) which is equivalent to a lumped
/// mass finite element (FE) scheme.
///
/////////////////////////////////////////////////////////////////////////////
class SGH : public Solver
{
public:

    char* MESH_FILE;

    SGH(); //SGH_Parameters& params, Solver* Solver_Pointer, std::shared_ptr<mesh_t> mesh_in, const int my_fea_module_index = 0);
    ~SGH();

    // initialize data for boundaries of the model and storage for boundary conditions and applied loads
    // void sgh_interface_setup(node_t& node, elem_t& elem, corner_t& corner){};

    void setup();

    void write_outputs(
        const mesh_t&              mesh,
        DViewCArrayKokkos<double>& node_coords,
        DViewCArrayKokkos<double>& node_vel,
        DViewCArrayKokkos<double>& node_mass,
        DViewCArrayKokkos<double>& elem_den,
        DViewCArrayKokkos<double>& elem_pres,
        DViewCArrayKokkos<double>& elem_stress,
        DViewCArrayKokkos<double>& elem_sspd,
        DViewCArrayKokkos<double>& elem_sie,
        DViewCArrayKokkos<double>& elem_vol,
        DViewCArrayKokkos<double>& elem_mass,
        DViewCArrayKokkos<size_t>& elem_mat_id,
        CArray<double>&            graphics_times,
        size_t&                    graphics_id,
        const double               time_value);

    // void cleanup_material_models();

    // void module_cleanup();

    void solve(CArrayKokkos<material_t>&  material,
               CArrayKokkos<boundary_t>&  boundary,
               mesh_t&                    mesh,
               DViewCArrayKokkos<double>& node_coords,
               DViewCArrayKokkos<double>& node_vel,
               DViewCArrayKokkos<double>& node_mass,
               DViewCArrayKokkos<double>& elem_den,
               DViewCArrayKokkos<double>& elem_pres,
               DViewCArrayKokkos<double>& elem_stress,
               DViewCArrayKokkos<double>& elem_sspd,
               DViewCArrayKokkos<double>& elem_sie,
               DViewCArrayKokkos<double>& elem_vol,
               DViewCArrayKokkos<double>& elem_div,
               DViewCArrayKokkos<double>& elem_mass,
               DViewCArrayKokkos<size_t>& elem_mat_id,
               DViewCArrayKokkos<double>& elem_statev,
               DViewCArrayKokkos<double>& corner_force,
               DViewCArrayKokkos<double>& corner_mass,
               double&                    time_value,
               const double               time_final,
               const double               dt_max,
               const double               dt_min,
               const double               dt_cfl,
               double&                    graphics_time,
               size_t                     graphics_cyc_ival,
               double                     graphics_dt_ival,
               const size_t               cycle_stop,
               const size_t               rk_num_stages,
               double                     dt,
               const double               fuzz,
               const double               tiny,
               const double               small,
               CArray<double>&            graphics_times,
               size_t&                    graphics_id);



    // **** Functions defined in boundary.cpp **** //
    void boundary_velocity(
        const mesh_t& mesh,
        const CArrayKokkos<boundary_t>& boundary,
        DViewCArrayKokkos<double>& node_vel,
        const double time_value);


    // **** Functions defined in energy_sgh.cpp **** //
    void update_energy(
        double rk_alpha,
        double dt,
        const mesh_t& mesh,
        const DViewCArrayKokkos<double>& node_vel,
        const DViewCArrayKokkos<double>& node_coords,
        DViewCArrayKokkos<double>&       elem_sie,
        const DViewCArrayKokkos<double>& elem_mass,
        const DViewCArrayKokkos<double>& corner_force);



    // **** Functions defined in eos.cpp **** //
    // NOTE: This should be moved up so multiple solvers can use the ideal gas equation
    KOKKOS_FUNCTION
    void ideal_gas(
        const DViewCArrayKokkos<double>& elem_pres,
        const DViewCArrayKokkos<double>& elem_stress,
        const size_t                     elem_gid,
        const size_t                     mat_id,
        const DViewCArrayKokkos<double>& elem_state_vars,
        const DViewCArrayKokkos<double>& elem_sspd,
        const double                     den,
        const double                     sie);

    // **** Functions defined in force_sgh.cpp **** //
    void get_force(
        const CArrayKokkos<material_t>&  material,
        const mesh_t&                    mesh,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        const DViewCArrayKokkos<double>& elem_den,
        const DViewCArrayKokkos<double>& elem_sie,
        const DViewCArrayKokkos<double>& elem_pres,
        const DViewCArrayKokkos<double>& elem_stress,
        const DViewCArrayKokkos<double>& elem_sspd,
        const DViewCArrayKokkos<double>& elem_vol,
        const DViewCArrayKokkos<double>& elem_div,
        const DViewCArrayKokkos<size_t>& elem_mat_id,
        DViewCArrayKokkos<double>&       corner_force,
        const double                     fuzz,
        const double                     small,
        const DViewCArrayKokkos<double>& elem_statev,
        const double                     dt,
        const double                     rk_alpha);


    void get_force_2D(
        const CArrayKokkos<material_t>&  material,
        const mesh_t&                    mesh,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        const DViewCArrayKokkos<double>& elem_den,
        const DViewCArrayKokkos<double>& elem_sie,
        const DViewCArrayKokkos<double>& elem_pres,
        const DViewCArrayKokkos<double>& elem_stress,
        const DViewCArrayKokkos<double>& elem_sspd,
        const DViewCArrayKokkos<double>& elem_vol,
        const DViewCArrayKokkos<double>& elem_div,
        const DViewCArrayKokkos<size_t>& elem_mat_id,
        DViewCArrayKokkos<double>&       corner_force,
        const double                     fuzz,
        const double                     small,
        const DViewCArrayKokkos<double>& elem_statev,
        const double                     dt,
        const double                     rk_alpha);



    // **** Functions defined in geometry.cpp **** //
    void update_position(
        double rk_alpha,
        double dt,
        const size_t num_dims,
        const size_t num_nodes,
        DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel);

    // NOTE: Consider pulling this up as well
    KOKKOS_FUNCTION
    void get_bmatrix(
        const ViewCArrayKokkos<double>&  B_matrix,
        const size_t                     elem_gid,
        const DViewCArrayKokkos<double>& node_coords,
        const ViewCArrayKokkos<size_t>&  elem_node_gids);

    // NOTE: Consider pulling this up as well
    void get_vol(
        const DViewCArrayKokkos<double>& elem_vol,
        const DViewCArrayKokkos<double>& node_coords,
        const mesh_t&                    mesh);


    // Exact volume for a hex element
    // NOTE: Consider pulling this up as well
    KOKKOS_FUNCTION
    void get_vol_hex(
        const DViewCArrayKokkos<double>& elem_vol,
        const size_t                     elem_gid,
        const DViewCArrayKokkos<double>& node_coords,
        const ViewCArrayKokkos<size_t>&  elem_node_gids);

    // true volume of a quad in RZ coords
    // NOTE: Consider pulling this up as well
    KOKKOS_FUNCTION
    void get_vol_quad(
        const DViewCArrayKokkos<double>& elem_vol,
        const size_t                     elem_gid,
        const DViewCArrayKokkos<double>& node_coords,
        const ViewCArrayKokkos<size_t>&  elem_node_gids);

    // NOTE: Consider pulling this up as well
    KOKKOS_FUNCTION
    void get_bmatrix2D(
        const ViewCArrayKokkos<double>&  B_matrix,
        const size_t                     elem_gid,
        const DViewCArrayKokkos<double>& node_coords,
        const ViewCArrayKokkos<size_t>&  elem_node_gids);


    // element facial area
    // NOTE: Consider pulling this up as well
    KOKKOS_FUNCTION
    double get_area_quad(
        const size_t elem_gid,
        const DViewCArrayKokkos<double>& node_coords,
        const ViewCArrayKokkos<size_t>&  elem_node_gids);

    // NOTE: Consider pulling this up as well
    KOKKOS_FUNCTION
    double heron(
        const double x1,
        const double y1,
        const double x2,
        const double y2,
        const double x3,
        const double y3);

    KOKKOS_FUNCTION
    void get_area_weights2D(
        const ViewCArrayKokkos<double>&  corner_areas,
        const size_t                     elem_gid,
        const DViewCArrayKokkos<double>& node_coords,
        const ViewCArrayKokkos<size_t>&  elem_node_gids);



    // **** Functions defined in momentum.cpp **** //
    void update_velocity(
        double rk_alpha,
        double dt,
        const mesh_t& mesh,
        DViewCArrayKokkos<double>&       node_vel,
        const DViewCArrayKokkos<double>& node_mass,
        const DViewCArrayKokkos<double>& corner_force);



    KOKKOS_FUNCTION
    void get_velgrad(
        ViewCArrayKokkos<double>&        vel_grad,
        const ViewCArrayKokkos<size_t>&  elem_node_gids,
        const DViewCArrayKokkos<double>& node_vel,
        const ViewCArrayKokkos<double>&  b_matrix,
        const double                     elem_vol,
        const size_t                     elem_gid);

    KOKKOS_FUNCTION
    void get_velgrad2D(
        ViewCArrayKokkos<double>&        vel_grad,
        const ViewCArrayKokkos<size_t>&  elem_node_gids,
        const DViewCArrayKokkos<double>& node_vel,
        const ViewCArrayKokkos<double>&  b_matrix,
        const double                     elem_vol,
        const double                     elem_area,
        const size_t                     elem_gid);


    void get_divergence(
        DViewCArrayKokkos<double>&       elem_div,
        const mesh_t                     mesh,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        const DViewCArrayKokkos<double>& elem_vol);

    void get_divergence2D(
        DViewCArrayKokkos<double>&       elem_div,
        const mesh_t                     mesh,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        const DViewCArrayKokkos<double>& elem_vol);

    KOKKOS_FUNCTION
    void decompose_vel_grad(
        ViewCArrayKokkos<double>&        D_tensor,
        ViewCArrayKokkos<double>&        W_tensor,
        const ViewCArrayKokkos<double>&  vel_grad,
        const ViewCArrayKokkos<size_t>&  elem_node_gids,
        const size_t                     elem_gid,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        const double                     vol);

    // **** Functions defined in properties.cpp **** //
    void update_state(
        const CArrayKokkos<material_t>&  material,
        const mesh_t&                    mesh,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        DViewCArrayKokkos<double>&       elem_den,
        DViewCArrayKokkos<double>&       elem_pres,
        DViewCArrayKokkos<double>&       elem_stress,
        DViewCArrayKokkos<double>&       elem_sspd,
        const DViewCArrayKokkos<double>& elem_sie,
        const DViewCArrayKokkos<double>& elem_vol,
        const DViewCArrayKokkos<double>& elem_mass,
        const DViewCArrayKokkos<size_t>& elem_mat_id,
        const DViewCArrayKokkos<double>& elem_statev,
        const double                     dt,
        const double                     rk_alpha);

    void update_state2D(
        const CArrayKokkos<material_t>&  material,
        const mesh_t&                    mesh,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        DViewCArrayKokkos<double>&       elem_den,
        DViewCArrayKokkos<double>&       elem_pres,
        DViewCArrayKokkos<double>&       elem_stress,
        DViewCArrayKokkos<double>&       elem_sspd,
        const DViewCArrayKokkos<double>& elem_sie,
        const DViewCArrayKokkos<double>& elem_vol,
        const DViewCArrayKokkos<double>& elem_mass,
        const DViewCArrayKokkos<size_t>& elem_mat_id,
        const DViewCArrayKokkos<double>& elem_statev,
        const double                     dt,
        const double                     rk_alpha);


    // **** Functions defined in time_integration.cpp **** //
    // NOTE: Consider pulling up
    void rk_init(
        DViewCArrayKokkos<double>& node_coords,
        DViewCArrayKokkos<double>& node_vel,
        DViewCArrayKokkos<double>& elem_sie,
        DViewCArrayKokkos<double>& elem_stress,
        const size_t               num_dims,
        const size_t               num_elems,
        const size_t               num_nodes);


    void get_timestep(
        mesh_t&                    mesh,
        DViewCArrayKokkos<double>& node_coords,
        DViewCArrayKokkos<double>& node_vel,
        DViewCArrayKokkos<double>& elem_sspd,
        DViewCArrayKokkos<double>& elem_vol,
        double                     time_value,
        const double               graphics_time,
        const double               time_final,
        const double               dt_max,
        const double               dt_min,
        const double               dt_cfl,
        double&                    dt,
        const double               fuzz);


    void get_timestep2D(
        mesh_t&                    mesh,
        DViewCArrayKokkos<double>& node_coords,
        DViewCArrayKokkos<double>& node_vel,
        DViewCArrayKokkos<double>& elem_sspd,
        DViewCArrayKokkos<double>& elem_vol,
        double                     time_value,
        const double               graphics_time,
        const double               time_final,
        const double               dt_max,
        const double               dt_min,
        const double               dt_cfl,
        double&                    dt,
        const double               fuzz);


    // **** Functions defined in user_mat.cpp **** //
    // NOTE: Pull up into high level
    KOKKOS_FUNCTION
    void user_eos_model(
        const DViewCArrayKokkos<double>& elem_pres,
        const DViewCArrayKokkos<double>& elem_stress,
        const size_t                     elem_gid,
        const size_t                     mat_id,
        const DViewCArrayKokkos<double>& elem_state_vars,
        const DViewCArrayKokkos<double>& elem_sspd,
        const double                     den,
        const double                     sie);


    KOKKOS_FUNCTION
    void user_strength_model(
        const DViewCArrayKokkos<double>& elem_pres,
        const DViewCArrayKokkos<double>& elem_stress,
        const size_t                     elem_gid,
        const size_t                     mat_id,
        const DViewCArrayKokkos<double>& elem_state_vars,
        const DViewCArrayKokkos<double>& elem_sspd,
        const double                     den,
        const double                     sie,
        const ViewCArrayKokkos<double>&  vel_grad,
        const ViewCArrayKokkos<size_t>&  elem_node_gids,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        const double                     vol,
        const double                     dt,
        const double                     rk_alpha);


    KOKKOS_FUNCTION
    void user_strength_model_vpsc(
        const DViewCArrayKokkos<double>& elem_pres,
        const DViewCArrayKokkos<double>& elem_stress,
        const size_t                     elem_gid,
        const size_t                     mat_id,
        const DViewCArrayKokkos<double>& elem_state_vars,
        const DViewCArrayKokkos<double>& elem_sspd,
        const double                     den,
        const double                     sie,
        const ViewCArrayKokkos<double>&  vel_grad,
        const ViewCArrayKokkos<size_t>&  elem_node_gids,
        const DViewCArrayKokkos<double>& node_coords,
        const DViewCArrayKokkos<double>& node_vel,
        const double                     vol,
        const double                     dt,
        const double                     rk_alpha);

    // **** Functions defined in user_mat_init.cpp **** //
    // NOTE: Pull up into high level
    void user_model_init(
        const DCArrayKokkos<double>& file_state_vars,
        const size_t                 num_state_vars,
        const size_t                 mat_id,
        const size_t                 num_elems);

};

#endif // end HEADER_H
