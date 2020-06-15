
#pragma once

#define PI          3.141592653589793238462
#define half_pi     1.570796326794897
#define degtorad    0.0174532925199
#define radtodeg    57.295779513082

class IK
{
public:
    struct IK_DATA {
        double d1; // Ground to q1
        double a2; // q1 to q2
        double d4; // q2 to wrist
        double d6; // wrist to gripper
    };

    IK(IK_DATA data) : m_IK(data) {}

    void inverseKinematics(double x, double y, double z, double t[3][3], double angles[6], bool flip);
    void forwardKinematics(double angles[6], double jointPos[6][3]);

    static double toRad(double a);
    static double toDeg(double a);

private:
    IK_DATA m_IK;
};
