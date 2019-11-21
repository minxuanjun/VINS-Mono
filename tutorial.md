1. 常用的代码片段

IMU 的中值积分

```
    double t = imu_msg->header.stamp.toSec();
    if (init_imu)
    {
        latest_time = t;
        init_imu = 0;
        return;
    }
    double dt = t - latest_time;
    latest_time = t;

    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    Eigen::Vector3d linear_acceleration{dx, dy, dz};

    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Eigen::Vector3d angular_velocity{rx, ry, rz};
    //　ａ　＝　R_g_imu *(am - ba) + g
    Eigen::Vector3d un_acc_0 = tmp_Q * (acc_0 - tmp_Ba) - estimator.g;
    // w = wm - bg
    // notice　这里角度的增量使用的是中值积分
    Eigen::Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - tmp_Bg;
    // R_j_i+1 = R_j_i * R_i_i+11
    tmp_Q = tmp_Q * Utility::deltaQ(un_gyr * dt);
    // 加速度进行中值积分
    Eigen::Vector3d un_acc_1 = tmp_Q * (linear_acceleration - tmp_Ba) - estimator.g;
    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
    //　核心公式
    tmp_P = tmp_P + dt * tmp_V + 0.5 * dt * dt * un_acc;
    tmp_V = tmp_V + dt * un_acc;

　　 // 存储上一时刻的加速度值, 
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
```

2. Eigen 的向量赋值

```
 Eigen::Matrix<double, 7, 1> xyz_uv_velocity;
                xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
```

3. Eigen 向量的单位化

```
Eigen::Vector3d f = it_per_frame.point.normalized();
```

4. 特征点三角化, 使用SVD进行求解

```
// 每个特征点的主导帧都是第一次观测到的第一帧
// 而basalt 的主导帧是中间的的一帧,也许这样更好, 但是basalt 只使用两帧初始化,而vins-mono是使用多帧进行三角化
Eigen::Vector3d t1 = Ps[imu_j] + Rs[imu_j] * tic[0]; 
Eigen::Matrix3d R1 = Rs[imu_j] * ric[0];
Eigen::Vector3d t = R0.transpose() * (t1 - t0);
Eigen::Matrix3d R = R0.transpose() * R1; 

Eigen::Matrix<double, 3, 4> P;
P.leftCols<3>() = R.transpose();
P.rightCols<1>() = -R.transpose() * t;
Eigen::Vector3d f = it_per_frame.point.normalized();
svd_A.row(svd_idx++) = f[0] * P.row(2) - f[2] * P.row(0);
svd_A.row(svd_idx++) = f[1] * P.row(2) - f[2] * P.row(1);
if (imu_i == imu_j) // 这句话没什么用
    continue;
Eigen::Vector4d svd_V = Eigen::JacobiSVD<Eigen::MatrixXd>(svd_A, Eigen::ComputeThinV).matrixV().rightCols<1>();
double svd_method = svd_V[2] / svd_V[3];
//it_per_id->estimated_depth = -b / A;
//it_per_id->estimated_depth = svd_V[2] / svd_V[3];

it_per_id.estimated_depth = svd_method;
//it_per_id->estimated_depth = INIT_DEPTH;

if (it_per_id.estimated_depth < 0.1)
{
    it_per_id.estimated_depth = INIT_DEPTH;
}
```

5. Eigen 取矩阵的列元素, ``` leftCols()``` 去矩阵的左边的几列元素, ```middleCols()``` 取中间的几列元素

```
Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobian(jacobians[i], n, size);
jacobian.setZero();
jacobian.leftCols(local_size) = marginalization_info->linearized_jacobians.middleCols(idx, local_size);
```

6. ceres　协方差矩阵转信息矩阵的平方根

6.2 LLT分解法

```
Eigen::Matrix<double, 15, 15> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 15, 15>>(pre_integration->covariance.inverse()).matrixL().transpose();
```

6.3 使用自对称矩阵分解法

实对称矩阵的特征向量是正交的,任何实对称矩阵都相似与一个对角阵
A = U*Sigma*U^{-1} = U*SIgma*U^{T}

[特征值和特征向量分解](http://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html)

```
Eigen::MatrixXd Amm = 0.5 * (A.block(0, 0, m, m) + A.block(0, 0, m, m).transpose());　
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Amm);

//ROS_ASSERT_MSG(saes.eigenvalues().minCoeff() >= -1e-4, "min eigenvalue %f", saes.eigenvalues().minCoeff());

Eigen::MatrixXd Amm_inv = saes.eigenvectors() * Eigen::VectorXd((saes.eigenvalues().array() > eps).select(saes.eigenvalues().array().inverse(), 0)).asDiagonal() * saes.eigenvectors().transpose();
```

7. VINS-Mono 的边缘化过程中的求逆运算

$H_{mm}^{-1}$的解法

method1:好像是[自对称矩阵分解法](http://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html),对较小的特征值能进行更加精细的处理，数值稳定性更好(VINS-Mono)
```
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Amm);

//ROS_ASSERT_MSG(saes.eigenvalues().minCoeff() >= -1e-4, "min eigenvalue %f", saes.eigenvalues().minCoeff());

Eigen::MatrixXd Amm_inv = saes.eigenvectors() * Eigen::VectorXd((saes.eigenvalues().array() > eps).select(saes.eigenvalues().array().inverse(), 0)).asDiagonal() * saes.eigenvectors().transpose();

```

method2: [LU分解](http://eigen.tuxfamily.org/dox/classEigen_1_1FullPivLU.html) (Basalt)

```

H_mm_inv = abs_H.bottomRightCorner(marg_size, marg_size)
        .fullPivLu()
        .solve(Eigen::MatrixXd::Identity(marg_size, marg_size));
```
８．VINS-Mono 的边缘化和Basalt 对比(需要大家讨论)
VINS-Mono: 边缘化KF0时，会边缘化其位姿和其主导的3d点，**但是其主导的3d点还留着（更换其主导帧)**，而basalt之间将其主导的3d点直接扔掉了，
边缘化次新帧Ｆ１，只边缘化其对应的先验，其主导的feature直接扔掉，同是将其他帧主导的3d点在次新帧的观测直接扔掉(为了稀疏行)
　