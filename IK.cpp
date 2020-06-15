#include <unistd.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include<vector>
#include <chrono>

#include "IK.h"

#define x_comp  0
#define y_comp  1
#define z_comp  2

using namespace::std;

double epsilon = 0.0001;

void IK::inverseKinematics(double x, double y, double z, double t[3][3], double angles[6], bool flip)
{
    double xc = x - m_IK.d6*t[0][2];
    double yc = y - m_IK.d6*t[1][2];
    double zc = z - m_IK.d6*t[2][2];

    angles[0] = atan2l(yc,xc);

    double D = ( pow(xc,2) + pow(yc,2) + pow((zc-m_IK.d1),2) - pow(m_IK.a2,2) - pow(m_IK.d4,2) ) / (2*m_IK.a2*m_IK.d4);
    angles[2] = atan2l(-sqrt(1 - pow(D,2)), D);

    double k1 = m_IK.a2+m_IK.d4*cosl(angles[2]);
    double k2 =         m_IK.d4*sinl(angles[2]);
    angles[1] = atan2l( (zc-m_IK.d1), sqrt(pow(xc,2) + pow(yc,2)) ) - atan2l(k2,k1) ;

    angles[2] += PI/2;

    double q1 = angles[0];
    double q2 = angles[1];
    double q3 = angles[2];

    double r11,r12,r13,r21,r22,r23,r31,r32,r33;
    r11=t[0][0]; r12=t[0][1]; r13=t[0][2];
    r21=t[1][0]; r22=t[1][1]; r23=t[1][2];
    r31=t[2][0]; r32=t[2][1]; r33=t[2][2];

    double ax =  r13*cosl(q1)*cosl(q2 + q3) + r23*cosl(q2 + q3)*sinl(q1) + r33*sinl(q2 + q3);
    double ay = -r23*cosl(q1)      + r13*sinl(q1);
    double az = -r33*cosl(q2 + q3) + r13*cosl(q1)*sinl(q2 + q3) + r23*sinl(q1)*sinl(q2 + q3);
    double sz = -r32*cosl(q2 + q3) + r12*cosl(q1)*sinl(q2 + q3) + r22*sinl(q1)*sinl(q2 + q3);
    double nz = -r31*cosl(q2 + q3) + r11*cosl(q1)*sinl(q2 + q3) + r21*sinl(q1)*sinl(q2 + q3);

    if (ax == 0) ax = epsilon;
    if (ay == 0) ay = epsilon;

	if (flip) {
		angles[3] = atan2(-ay,-ax);
		angles[4] = atan2(-sqrt(ax*ax+ay*ay),az);
		angles[5] = atan2(-sz,nz);
	}
	else {
		angles[3] = atan2(ay,ax);
		angles[4] = atan2(sqrt(ax*ax+ay*ay),az);
		angles[5] = atan2(sz,-nz);
	}
}

void IK::forwardKinematics(double angles[6], double jointPos[6][3])
{
    double q1,q2,q3,q4,q5;
    q1=angles[0]; q2=angles[1]; q3=angles[2]; q4=angles[3]; q5=angles[4];

    jointPos[0][x_comp] = 0;
    jointPos[0][y_comp] = 0;
    jointPos[0][z_comp] = m_IK.d1;

    jointPos[1][x_comp]  = m_IK.a2*cos(q1)*cos(q2);
    jointPos[1][y_comp]  = m_IK.a2*cos(q2)*sin(q1);
    jointPos[1][z_comp]  = m_IK.d1+m_IK.a2*sin(q2);

    jointPos[2][x_comp] = 0;
    jointPos[2][y_comp] = 0;
    jointPos[2][z_comp] = 0;
    // jointPos[2][x_comp] = cos(q1)*(m_IK.a2*cos(q2) + (m_IK.d4/2.0)*sin(q2+q3));
    // jointPos[2][y_comp] = sin(q1)*(m_IK.a2*cos(q2) + (m_IK.d4/2.0)*sin(q2+q3));
    // jointPos[2][z_comp] = m_IK.d1 - (m_IK.d4/2.0)*cos(q2+q3)+m_IK.a2*sin(q2);

    jointPos[3][x_comp] = cos(q1)*(m_IK.a2*cos(q2) + m_IK.d4*sin(q2+q3));
    jointPos[3][y_comp] = sin(q1)*(m_IK.a2*cos(q2) + m_IK.d4*sin(q2+q3));
    jointPos[3][z_comp] = m_IK.d1 - m_IK.d4*cos(q2+q3) + m_IK.a2*sin(q2);

    jointPos[4][x_comp] = 0;
    jointPos[4][y_comp] = 0;
    jointPos[4][z_comp] = 0;
    // jointPos[4][x_comp] = (m_IK.d6/2.0)*sin(q1)*sin(q4)*sin(q5) + cos(q1)*(m_IK.a2*cos(q2) + (m_IK.d4 + (m_IK.d6/2.0)*cos(q5))*sin(q2 + q3) + (m_IK.d6/2.0)*cos(q2 + q3)*cos(q4)*sin(q5));
    // jointPos[4][y_comp] = cos(q3)*(m_IK.d4 + (m_IK.d6/2.0)*cos(q5))*sin(q1)*sin(q2) - (m_IK.d6/2.0)*(cos(q4)*sin(q1)*sin(q2)*sin(q3) + cos(q1)*sin(q4))*sin(q5) + cos(q2)*sin(q1)*(m_IK.a2 + (m_IK.d4 + (m_IK.d6/2.0)*cos(q5))*sin(q3) + (m_IK.d6/2.0)*cos(q3)*cos(q4)*sin(q5));
    // jointPos[4][z_comp] = m_IK.d1 - cos(q2 + q3)*(m_IK.d4 + (m_IK.d6/2.0)*cos(q5)) + m_IK.a2*sin(q2) + (m_IK.d6/2.0)*cos(q4)*sin(q2 + q3)*sin(q5);

    jointPos[5][x_comp] = m_IK.d6*sin(q1)*sin(q4)*sin(q5) + cos(q1)*(m_IK.a2*cos(q2) + (m_IK.d4 + m_IK.d6*cos(q5))*sin(q2 + q3) + m_IK.d6*cos(q2 + q3)*cos(q4)*sin(q5));
    jointPos[5][y_comp] = cos(q3)*(m_IK.d4 + m_IK.d6*cos(q5))*sin(q1)*sin(q2) - m_IK.d6*(cos(q4)*sin(q1)*sin(q2)*sin(q3) + cos(q1)*sin(q4))*sin(q5) + cos(q2)*sin(q1)*(m_IK.a2 + (m_IK.d4 + m_IK.d6*cos(q5))*sin(q3) + m_IK.d6*cos(q3)*cos(q4)*sin(q5));
    jointPos[5][z_comp] = m_IK.d1 - cos(q2 + q3)*(m_IK.d4 + m_IK.d6*cos(q5)) + m_IK.a2*sin(q2) + m_IK.d6*cos(q4)*sin(q2 + q3)*sin(q5);
}

double IK::toRad(double a) { return a * degtorad; }
double IK::toDeg(double a) { return a * radtodeg; }

void Dot(double m1[3][3], double m2[3][3], double m3[3][3])
{
    for (int y = 0; y < 3; ++y) {
        for (int x = 0; x < 3; ++x) {
            double sum = 0;
            for (int i = 0; i < 3; ++i)
                sum += m1[y][i] * m2[i][x];

            m3[y][x] = sum;
        }
    }
}
