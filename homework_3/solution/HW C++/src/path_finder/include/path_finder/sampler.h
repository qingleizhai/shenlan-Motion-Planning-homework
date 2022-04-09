/*
Copyright (C) 2022 Hongkai Ye (kyle_yeh@163.com)
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
#ifndef _BIAS_SAMPLER_
#define _BIAS_SAMPLER_

#include <ros/ros.h>
#include <Eigen/Eigen>
#include <random>

class BiasSampler
{
public:
  BiasSampler()
  {
    std::random_device rd;
    gen_ = std::mt19937_64(rd());
    uniform_rand_ = std::uniform_real_distribution<double>(0.0, 1.0);
    normal_rand_ = std::normal_distribution<double>(0.0, 1.0);
    range_.setZero();
    origin_.setZero();
  };

  void setSamplingRange(const Eigen::Vector3d origin, const Eigen::Vector3d range)
  {
    origin_ = origin;
    range_ = range;
  }

  void samplingOnce(Eigen::Vector3d &sample)
  {
    sample[0] = uniform_rand_(gen_);
    sample[1] = uniform_rand_(gen_);
    sample[2] = uniform_rand_(gen_);
    sample.array() *= range_.array();
    sample += origin_;
  };

  void samplingOnceInEllipse(Eigen::Vector3d &sample, double c_max, const Eigen::Vector3d &s, const Eigen::Vector3d &g)
  {

    Eigen::Vector3d x_center = (s + g) / 2;

    double c_min = (s - g).norm();

    Eigen::Matrix3d L;
    calcDiagMatrix(L, c_max, c_min);
    sample = L*samplingUnitNBall();

    
    Eigen::Matrix3d C;
    calcRotationToWorldFrame(C, s, g);

    Eigen::Vector3d sample1 = C * sample;

    sample = sample1 + x_center ;
  }

  void calcDiagMatrix(Eigen::Matrix3d & diag, double c_best, double c_min)
  {
    diag.setZero();
    diag << c_best / 2, 0, 0,
            0,sqrt(c_best * c_best - c_min * c_min) / 2.0 ,0,
            0, 0, sqrt(c_best * c_best - c_min * c_min) / 2.0;
  }

  void calcRotationToWorldFrame(Eigen::Matrix3d & C, const Eigen::Vector3d &s, const Eigen::Vector3d &g)
  {
    auto dir = (g - s).normalized();

    C.col(0) = dir;
    dir[2] = 0.0;
    C.col(1) = Eigen::AngleAxisd(0.5 * M_PI, Eigen::Vector3d::UnitZ()) * dir.normalized(); // project to the x-y plane and then rotate 90 degree;
    C.col(2) = C.col(0).cross(C.col(1));

    //double beta = std::atan()


  }

  Eigen::Vector3d samplingUnitNBall()
  {
    Eigen::Vector3d sample;
    
    sample[0] = normal_rand_(gen_);
    sample[1] = normal_rand_(gen_);
    sample[2] = normal_rand_(gen_);
    double r = pow(uniform_rand_(gen_), 0.33333);
    
    return r * sample.normalized();
    
  }

  // (0.0 - 1.0)
  double getUniRandNum()
  {
    return uniform_rand_(gen_);
  }

private:
  Eigen::Vector3d range_, origin_;
  std::mt19937_64 gen_;
  std::uniform_real_distribution<double> uniform_rand_;
  std::normal_distribution<double> normal_rand_;
};

#endif
