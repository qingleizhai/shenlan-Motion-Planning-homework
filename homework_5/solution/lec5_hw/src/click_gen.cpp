#include "lec5_hw/visualizer.hpp"
#include "lec5_hw/trajectory.hpp"

#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <vector>

struct Config
{
    std::string targetTopic;
    double clickHeight;
    std::vector<double> initialVel;
    std::vector<double> initialAcc;
    std::vector<double> terminalVel;
    std::vector<double> terminalAcc;
    double allocationSpeed;
    double allocationAcc;
    int maxPieceNum;

    Config(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("ClickHeight", clickHeight);
        nh_priv.getParam("InitialVel", initialVel);
        nh_priv.getParam("InitialAcc", initialAcc);
        nh_priv.getParam("TerminalVel", terminalVel);
        nh_priv.getParam("TerminalAcc", terminalAcc);
        nh_priv.getParam("AllocationSpeed", allocationSpeed);
        nh_priv.getParam("AllocationAcc", allocationAcc);
        nh_priv.getParam("MaxPieceNum", maxPieceNum);
    }
};

double timeTrapzVel(const double dist,
                    const double vel,
                    const double acc)
{
    const double t = vel / acc;
    const double d = 0.5 * acc * t * t;

    if (dist < d + d)
    {
        return 2.0 * sqrt(dist / acc);
    }
    else
    {
        return 2.0 * t + (dist - 2.0 * d) / vel;
    }
}

void minimumJerkTrajGen(
    // Inputs:
    const int pieceNum,
    const Eigen::Vector3d &initialPos,
    const Eigen::Vector3d &initialVel,
    const Eigen::Vector3d &initialAcc,
    const Eigen::Vector3d &terminalPos,
    const Eigen::Vector3d &terminalVel,
    const Eigen::Vector3d &terminalAcc,
    const Eigen::Matrix3Xd &intermediatePositions,
    const Eigen::VectorXd &timeAllocationVector,
    // Outputs:
    Eigen::MatrixX3d &coefficientMatrix)
{
    // coefficientMatrix is a matrix with 6*piece num rows and 3 columes
    // As for a polynomial c0+c1*t+c2*t^2+c3*t^3+c4*t^4+c5*t^5,
    // each 6*3 sub-block of coefficientMatrix is
    // --              --
    // | c0_x c0_y c0_z |
    // | c1_x c1_y c1_z |
    // | c2_x c2_y c2_z |
    // | c3_x c3_y c3_z |
    // | c4_x c4_y c4_z |
    // | c5_x c5_y c5_z |
    // --              --
    // Please computed coefficientMatrix of the minimum-jerk trajectory
    // in this function

    // ------------------------ Put your solution below ------------------------

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(6 * pieceNum, 6 * pieceNum);
    Eigen::MatrixX3d B = Eigen::MatrixXd::Zero(6 * pieceNum, 3);

    // initial p,v,a constraint
    M(0, 0) = 1.0;
    M(1, 1) = 1.0;
    M(2, 2) = 2.0;
    B.row(0) = initialPos.transpose();
    B.row(1) = initialVel.transpose();
    B.row(2) = initialAcc.transpose();


    for (int i = 0; i < pieceNum - 1; i++) {
        
        B.row(6 * i + 3) = intermediatePositions.col(i).transpose();

        // position constaint  
        M(6 * i + 3, 6 * i) = 1.0;
        M(6 * i + 3, 6 * i + 1) = timeAllocationVector(i);
        M(6 * i + 3, 6 * i + 2) = std::pow(timeAllocationVector(i), 2);
        M(6 * i + 3, 6 * i + 3) = std::pow(timeAllocationVector(i), 3);
        M(6 * i + 3, 6 * i + 4) = std::pow(timeAllocationVector(i), 4);
        M(6 * i + 3, 6 * i + 5) = std::pow(timeAllocationVector(i), 5);

        // position continuity constraint
        M(6 * i + 4, 6 * i) = 1.0;
        M(6 * i + 4, 6 * i + 1) = timeAllocationVector(i);
        M(6 * i + 4, 6 * i + 2) = std::pow(timeAllocationVector(i), 2);
        M(6 * i + 4, 6 * i + 3) = std::pow(timeAllocationVector(i), 3);
        M(6 * i + 4, 6 * i + 4) = std::pow(timeAllocationVector(i), 4);
        M(6 * i + 4, 6 * i + 5) = std::pow(timeAllocationVector(i), 5);
        M(6 * i + 4, 6 * i + 6) = -1.0;

        // velocity continuity constraint
        M(6 * i + 5, 6 * i + 1) = 1.0;
        M(6 * i + 5, 6 * i + 2) = 2 * timeAllocationVector(i);
        M(6 * i + 5, 6 * i + 3) = 3 * std::pow(timeAllocationVector(i), 2);
        M(6 * i + 5, 6 * i + 4) = 4 * std::pow(timeAllocationVector(i), 3);
        M(6 * i + 5, 6 * i + 5) = 5 * std::pow(timeAllocationVector(i), 4);
        M(6 * i + 5, 6 * i + 7) = -1.0;

        // acceleration continuity constraint
        M(6 * i + 6, 6 * i + 2) = 2.0;
        M(6 * i + 6, 6 * i + 3) = 6 * timeAllocationVector(i);
        M(6 * i + 6, 6 * i + 4) = 12 * std::pow(timeAllocationVector(i), 2);
        M(6 * i + 6, 6 * i + 5) = 20 * std::pow(timeAllocationVector(i), 3);
        M(6 * i + 6, 6 * i + 8) = -2.0;

         // jerk continuity constraint
        M(6 * i + 7, 6 * i + 3) = 6.0;
        M(6 * i + 7, 6 * i + 4) = 24.0 * timeAllocationVector(i);
        M(6 * i + 7, 6 * i + 5) = 60.0 * std::pow(timeAllocationVector(i), 2);
        M(6 * i + 7, 6 * i + 9) = -6.0;

        // snap continuity constraint
        M(6 * i + 8, 6 * i + 4) = 24.0;
        M(6 * i + 8, 6 * i + 5) = 120.0 * timeAllocationVector(i);
        M(6 * i + 8, 6 * i + 10) = -24.0;
    }

    M(6 * pieceNum - 3, 6 * pieceNum - 6) = 1.0;
    M(6 * pieceNum - 3, 6 * pieceNum - 5) = timeAllocationVector(pieceNum - 1);
    M(6 * pieceNum - 3, 6 * pieceNum - 4) = std::pow(timeAllocationVector(pieceNum - 1), 2);
    M(6 * pieceNum - 3, 6 * pieceNum - 3) = std::pow(timeAllocationVector(pieceNum - 1), 3);
    M(6 * pieceNum - 3, 6 * pieceNum - 2) = std::pow(timeAllocationVector(pieceNum - 1), 4);
    M(6 * pieceNum - 3, 6 * pieceNum - 1) = std::pow(timeAllocationVector(pieceNum - 1), 5);
    M(6 * pieceNum - 2, 6 * pieceNum - 5) = 1.0;
    M(6 * pieceNum - 2, 6 * pieceNum - 4) = 2 * timeAllocationVector(pieceNum - 1);
    M(6 * pieceNum - 2, 6 * pieceNum - 3) = 3 * std::pow(timeAllocationVector(pieceNum - 1), 2);
    M(6 * pieceNum - 2, 6 * pieceNum - 2) = 4 * std::pow(timeAllocationVector(pieceNum - 1), 3);
    M(6 * pieceNum - 2, 6 * pieceNum - 1) = 5 * std::pow(timeAllocationVector(pieceNum - 1), 4);
    M(6 * pieceNum - 1, 6 * pieceNum - 4) = 2;
    M(6 * pieceNum - 1, 6 * pieceNum - 3) = 6 * timeAllocationVector(pieceNum - 1);
    M(6 * pieceNum - 1, 6 * pieceNum - 2) = 12 * std::pow(timeAllocationVector(pieceNum - 1), 2);
    M(6 * pieceNum - 1, 6 * pieceNum - 1) = 20 * std::pow(timeAllocationVector(pieceNum - 1), 3);

    B.row(6 * pieceNum - 3) = terminalPos.transpose();
    B.row(6 * pieceNum - 2) = terminalVel.transpose();
    B.row(6 * pieceNum - 1) = terminalAcc.transpose();

    
    coefficientMatrix =  M.inverse() * B;
    // ------------------------ Put your solution above ------------------------
}

class ClickGen
{
private:
    Config config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;

    Eigen::Matrix3Xd positions;
    Eigen::VectorXd times;
    int positionNum;
    Trajectory<5> traj;

public:
    ClickGen(const Config &conf,
             ros::NodeHandle &nh_)
        : config(conf),
          nh(nh_),
          visualizer(nh),
          positions(3, config.maxPieceNum + 1),
          times(config.maxPieceNum),
          positionNum(0)
    {
        targetSub = nh.subscribe(config.targetTopic, 1,
                                 &ClickGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        if (positionNum > config.maxPieceNum)
        {
            positionNum = 0;
            traj.clear();
        }

        positions(0, positionNum) = msg->pose.position.x;
        positions(1, positionNum) = msg->pose.position.y;
        positions(2, positionNum) = std::fabs(msg->pose.orientation.z) * config.clickHeight;

        if (positionNum > 0)
        {
            const double dist = (positions.col(positionNum) - positions.col(positionNum - 1)).norm();
            times(positionNum - 1) = timeTrapzVel(dist, config.allocationSpeed, config.allocationAcc);
        }

        ++positionNum;

        if (positionNum > 1)
        {
            const int pieceNum = positionNum - 1;
            const Eigen::Vector3d initialPos = positions.col(0);
            const Eigen::Vector3d initialVel(config.initialVel[0], config.initialVel[1], config.initialVel[2]);
            const Eigen::Vector3d initialAcc(config.initialAcc[0], config.initialAcc[1], config.initialAcc[2]);
            const Eigen::Vector3d terminalPos = positions.col(pieceNum);
            const Eigen::Vector3d terminalVel(config.terminalVel[0], config.terminalVel[1], config.terminalVel[2]);
            const Eigen::Vector3d terminalAcc(config.terminalAcc[0], config.terminalAcc[1], config.terminalAcc[2]);
            const Eigen::Matrix3Xd intermediatePositions = positions.middleCols(1, pieceNum - 1);
            const Eigen::VectorXd timeAllocationVector = times.head(pieceNum);

            Eigen::MatrixX3d coefficientMatrix = Eigen::MatrixXd::Zero(6 * pieceNum, 3);

            minimumJerkTrajGen(pieceNum,
                               initialPos, initialVel, initialAcc,
                               terminalPos, terminalVel, terminalAcc,
                               intermediatePositions,
                               timeAllocationVector,
                               coefficientMatrix);

            traj.clear();
            traj.reserve(pieceNum);
            for (int i = 0; i < pieceNum; i++)
            {
                traj.emplace_back(timeAllocationVector(i),
                                  coefficientMatrix.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
            }
             std::cout << coefficientMatrix << std::endl;
        }

        visualizer.visualize(traj, positions.leftCols(positionNum));

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "click_gen_node");
    ros::NodeHandle nh_;
    ClickGen clickGen(Config(ros::NodeHandle("~")), nh_);
    ros::spin();
    return 0;
}
