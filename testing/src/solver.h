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

#ifndef FIERRO_SOLVER_H
#define FIERRO_SOLVER_H

// #include "utilities.h"
#include "matar.h"
// #include "elements.h"
// #include "node_combination.h"
// #include "Simulation_Parameters/Simulation_Parameters.h"
// #include "FEA_Module.h"


// #include "MeshBuilder.h"
#include <map>
#include <memory>

using namespace mtr;

// forward declarations
namespace swage
{
class mesh_t;
} // namespace swage


class Solver
{
public:

   

    Solver(char* MESH);//Simulation_Parameters& _simparam);
    
    virtual ~Solver();

    virtual void setup() {}

    virtual void run() = 0;

    void solver_setup() {}

    void solver_finalize() {}

    void exit_solver(int status);

    // debug and system functions/variables
    double CPU_Time();
    void init_clock();

    double initial_CPU_time, communication_time, dev2host_time, host2dev_time, output_time;

    // class Simulation_Parameters *simparam;
    // Simulation_Parameters simparam;

    // set of enabled FEA modules
    // std::vector<FEA_MODULE_TYPE> fea_module_types;
    // std::vector<FEA_Module*>     fea_modules;
    // std::vector<bool> fea_modules_modal_analysis;
    // std::set<FEA_MODULE_TYPE> fea_module_must_read;
    int nfea_modules;
    // int displacement_module;

    bool finalize_flag = false;

    // elements::Element2D* elem2D;
    // elements::Element3D* elem;


    
};

#endif // end Header Guard