#include <hw_tool.h>
#include <unsupported/Eigen/Polynomials>

using namespace std;
using namespace Eigen;

void Homeworktool::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
}

void Homeworktool::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      
    
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

bool Homeworktool::isObsFree(const double coord_x, const double coord_y, const double coord_z)
{
    Vector3d pt;
    Vector3i idx;
    
    pt(0) = coord_x;
    pt(1) = coord_y;
    pt(2) = coord_z;
    idx = coord2gridIndex(pt);

    int idx_x = idx(0);
    int idx_y = idx(1);
    int idx_z = idx(2);

    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

Vector3d Homeworktool::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i Homeworktool::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d Homeworktool::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

double Homeworktool::OptimalBVP(Eigen::Vector3d _start_position,Eigen::Vector3d _start_velocity,Eigen::Vector3d _target_position)
{
    double optimal_cost = 100000; // this just to initial the optimal_cost, you can delete it 
    /*
                    



    STEP 2: go to the hw_tool.cpp and finish the function Homeworktool::OptimalBVP
    the solving process has been given in the document

    because the final point of trajectory is the start point of OBVP, so we input the pos,vel to the OBVP

    after finish Homeworktool::OptimalBVP, the Trajctory_Cost will record the optimal cost of this trajectory


    */

    auto p_x_0 = _start_position.x();
    auto p_y_0 = _start_position.y();
    auto p_z_0 = _start_position.z();

    auto p_x_0_sq = p_x_0 * p_x_0;
    auto p_y_0_sq = p_y_0 * p_y_0;
    auto p_z_0_sq = p_z_0 * p_z_0;
    
    auto v_x_0 = _start_velocity.x();
    auto v_y_0 = _start_velocity.y();
    auto v_z_0 = _start_velocity.z();

    auto v_x_0_sq = v_x_0 * v_x_0;
    auto v_y_0_sq = v_y_0 * v_y_0;
    auto v_z_0_sq = v_z_0 * v_z_0;

    auto p_x_f = _target_position.x();
    auto p_y_f = _target_position.y();
    auto p_z_f = _target_position.z();

    auto p_x_f_sq = p_x_f * p_x_f;
    auto p_y_f_sq = p_y_f * p_y_f;
    auto p_z_f_sq = p_z_f * p_z_f;

    
    auto v_x_f = _start_velocity.x();
    auto v_y_f = _start_velocity.y();
    auto v_z_f = _start_velocity.z();

    auto v_x_f_sq = v_x_f * p_x_f;
    auto v_y_f_sq = v_y_f * p_y_f;
    auto v_z_f_sq = v_z_f * p_z_f;

    auto cal_cost = [&](const double T) {
        double cost = (
            std::pow(T,4) + T * T * (4*v_x_0_sq + 4*v_x_0*v_x_f + 4*v_x_f_sq + 4*v_y_0_sq + 4*v_y_0*v_y_f + 4*v_y_f_sq + 4*v_z_0_sq + 4*v_z_0*v_z_f + 4*v_z_f_sq) + 
            T * (12*p_x_f*v_x_0 + 12*p_x_0*v_x_f - 12*p_x_0*v_x_0 - 12*p_x_f*v_x_f + 12*p_y_0*v_y_0 + 12*p_y_0*v_y_f - 12*p_y_f*v_y_0 - 12*p_y_f*v_y_f + 12*p_z_0*v_z_0 + 12*p_z_0*v_z_f - 12*p_z_f*v_z_0 - 12*p_z_f*v_z_f) - 
            12*p_x_0_sq - 24*p_x_0*p_x_f + 12*p_x_f_sq + 12*p_y_0_sq - 24*p_y_0*p_y_f + 12*p_y_f_sq + 12*p_z_0_sq - 24*p_z_0*p_z_f + 12*p_z_f_sq

        );
        return cost / std::pow(T, 3);
    };

    /*J= (T^4 + 
    (4*v_x_0^2 + 4*v_x_0*v_x_f + 4*v_x_f^2 + 4*v_y_0^2 + 4*v_y_0*v_y_f + 4*v_y_f^2 + 4*v_z_0^2 + 4*v_z_0*v_z_f + 4*v_z_f^2)*T^2 + 
    (12*p_x_0*v_x_0 + 12*p_x_0*v_x_f - 12*p_x_f*v_x_0 - 12*p_x_f*v_x_f + 12*p_y_0*v_y_0 + 12*p_y_0*v_y_f - 12*p_y_f*v_y_0 - 12*p_y_f*v_y_f + 12*p_z_0*v_z_0 + 12*p_z_0*v_z_f - 12*p_z_f*v_z_0 - 12*p_z_f*v_z_f)*T + 
    12*p_x_0^2 - 24*p_x_0*p_x_f + 12*p_x_f^2 + 12*p_y_0^2 - 24*p_y_0*p_y_f + 12*p_y_f^2 + 12*p_z_0^2 - 24*p_z_0*p_z_f + 12*p_z_f^2)/T^3
    */

   // derivative = 0
   /*
       c = [1, 
       - 4*v_x_0^2 - 4*v_x_0*v_x_f - 4*v_x_f^2 - 4*v_y_0^2 - 4*v_y_0*v_y_f - 4*v_y_f^2 - 4*v_z_0^2 - 4*v_z_0*v_z_f - 4*v_z_f^2, 
       24*p_x_f*v_x_0 - 24*p_x_0*v_x_f - 24*p_x_0*v_x_0 + 24*p_x_f*v_x_f - 24*p_y_0*v_y_0 - 24*p_y_0*v_y_f + 24*p_y_f*v_y_0 + 24*p_y_f*v_y_f - 24*p_z_0*v_z_0 - 24*p_z_0*v_z_f + 24*p_z_f*v_z_0 + 24*p_z_f*v_z_f, 
       - 36*p_x_0^2 + 72*p_x_0*p_x_f - 36*p_x_f^2 - 36*p_y_0^2 + 72*p_y_0*p_y_f - 36*p_y_f^2 - 36*p_z_0^2 + 72*p_z_0*p_z_f - 36*p_z_f^2]
 
      t =[T^6, T^4, T^3, T^2]
   */

    Eigen::VectorXd c(7);

    c << 0.0,
         0.0,
         -36*p_x_0_sq + 72*p_x_0*p_x_f - 36*p_x_f_sq - 36*p_y_0_sq + 72*p_y_0*p_y_f - 36*p_y_f_sq - 36*p_z_0_sq + 72*p_z_0*p_z_f - 36*p_z_f_sq,
         24*p_x_0*v_x_0 - 24*p_x_0*v_x_f - 24*p_x_f*v_x_0 + 24*p_x_f*v_x_f - 24*p_y_0*v_y_0 - 24*p_y_0*v_y_f + 24*p_y_f*v_y_0 + 24*p_y_f*v_y_f - 24*p_z_0*v_z_0 - 24*p_z_0*v_z_f + 24*p_z_f*v_z_0 + 24*p_z_f*v_z_f,
         -4*v_x_0_sq - 4*v_x_0*v_x_f - 4*v_x_f_sq - 4*v_y_0_sq - 4*v_y_0*v_y_f - 4*v_y_f_sq - 4*v_z_0_sq - 4*v_z_0*v_z_f - 4*v_z_f_sq,
         0.0,
         1.0;

    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    solver.compute(c);

    const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &rt = solver.roots();


    for (int i = 0; i < rt.rows(); ++ i) {
        if (rt(i).real() - 0.0 > 1e-6 && std::abs(rt(i).imag()) - 0.0 < 1e-6) {

            auto curr_cost = cal_cost(rt(i).real());
            optimal_cost = curr_cost < optimal_cost ? curr_cost : optimal_cost;

            std::cout << "curr_cost : " << curr_cost << std::endl;
        }

        
    }

    return optimal_cost;
}
