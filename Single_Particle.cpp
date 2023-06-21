#include <iostream>
#include <bits/stdc++.h>
#include <fstream>
using namespace std;

class quaternion // Defining the quaternion for ease of calculations
{
    public: float s, vx, vy, vz;

    quaternion()
    {
        s=0;
        vx=0;
        vy=0;
        vz=0;
    }

    quaternion(float a, float b, float c, float d)
    {
        s = a;
        vx = b;
        vy = c;
        vz = d;
    }

    public: quaternion product(quaternion q1, quaternion q2)
    {
        quaternion p;
        p.s = (q1.s)*(q2.s)-(q1.vx)*(q2.vx)-(q1.vy)*(q2.vy)-(q1.vz)*(q2.vz);
        p.vx = (q1.s)*(q2.vx)+(q2.s)*(q1.vx)+((q1.vy)*(q2.vz)-(q1.vz)*(q2.vy));
        p.vy = (q1.s)*(q2.vy)+(q2.s)*(q1.vy)+((q1.vz)*(q2.vx)-(q1.vx)*(q2.vz));
        p.vz = (q1.s)*(q2.vz)+(q2.s)*(q1.vz)+((q1.vx)*(q2.vy)-(q1.vy)*(q2.vx));
        return p;
    }

    public: quaternion T(quaternion q1)
    {
        quaternion p;
        p.s = q1.s;
        p.vx = -q1.vx;
        p.vy = -q1.vy;
        p.vz = -q1.vz;

        return p;
    }

    public: quaternion sum(quaternion q1, quaternion q2)
    {
        quaternion s;
        s.s = q1.s+q2.s;
        s.vx = q1.vx+q2.vx;
        s.vy = q1.vy+q2.vy;
        s.vz = q1.vz+q2.vz;

        return s;
    }

    public: quaternion diff(quaternion q1, quaternion q2)
    {
        quaternion d;
        d.s = q1.s - q2.s;
        d.vx = q1.vx - q2.vx;
        d.vy = q1.vy - q2.vy;
        d.vz = q1.vz - q2.vz;
        return d;
    }

    public: void display(quaternion q)
    {
        cout <<"("<<q.s<<","<<q.vx<<","<<q.vy<<","<<q.vz<<")"<<"\n";
        return;
    }
};


class sphere // Defining a Sphere
{
    public: quaternion rs, rb, F;
    public: float R, m;

    sphere()
    {
        R=0;
    }

    sphere(float x, float y, float z, float r, float mass)
    {
        rb.vx = x;
        rb.vy = y;
        rb.vz = z;
        R = r;
        m = mass;
    }
    
};

class particle // Defining a Particle
{

    public: quaternion rcm, q, wb, vcm, dq, T, F, acm, dw, I;
    public: int N_spheres=3;
    public: float m, k, eta;
    public: sphere A = sphere(0,1,0,0.9,1);
    public: sphere B = sphere(-0.866,-0.5,0,0.9,1);
    public: sphere C = sphere(0.866,-0.5,0,0.9,1);

    particle(float mass, float kn, float et, float x, float y, float z, float vx, float vy, float vz, float e1, float e2, float e3, float e4)
    {
        m = mass;
        k = kn;
        eta = et;
        rcm.vx = x;
        rcm.vy = y;
        rcm.vz = z;
        vcm.vx = vx;
        vcm.vy = vy;
        vcm.vz = vz;
        wb.vx = 0;
        wb.vy = 0;
        wb.vz = 0;

        I.vx = 3*2*A.m*A.R*A.R/5 + (A.rb.vy)*(A.rb.vy) + (A.rb.vz)*(A.rb.vz) + (B.rb.vy)*(B.rb.vy) + (B.rb.vz)*(B.rb.vz) + (C.rb.vy)*(C.rb.vy) + (C.rb.vz)*(C.rb.vz);
        I.vy = 3*2*A.m*A.R*A.R/5 + (A.rb.vx)*(A.rb.vx) + (A.rb.vy)*(A.rb.vy) + (B.rb.vx)*(B.rb.vx) + (B.rb.vy)*(B.rb.vy) + (C.rb.vx)*(C.rb.vx) + (C.rb.vy)*(C.rb.vy);
        I.vz = 3*2*A.m*A.R*A.R/5 + (A.rb.vz)*(A.rb.vz) + (A.rb.vx)*(A.rb.vx) + (B.rb.vz)*(B.rb.vz) + (B.rb.vx)*(B.rb.vx) + (C.rb.vz)*(C.rb.vz) + (C.rb.vx)*(C.rb.vx);

        q.s = e1;
        q.vx = e2;
        q.vy = e3;
        q.vz = e4;
        
    }

    public: void computeNetF()
    {
        quaternion NF;

        NF = NF.sum(NF, A.F);
        NF = NF.sum(NF, B.F);
        NF = NF.sum(NF, C.F);
        F  = NF;
        
    }


    public: void computeNetT()
    {
        quaternion NT;
        NT = NT.sum(NT, NT.product(A.rb, A.F));
        NT = NT.sum(NT, NT.product(B.rb, B.F));
        NT = NT.sum(NT, NT.product(C.rb, C.F));
        NT.s = 0;
        T = NT;
    }


};


class simDom // Defining the Simulation Domain
{

    public: float x1,y1,z1;
    public: float x2,y2,z2;
    simDom(float a,float b,float c,float d,float e,float f)
    {
        x1 = a;
        y1 = b;
        z1 = c;
        x2 = d;
        y2 = e;
        z2 = f;
    }

};


float calculateForce(float c, particle Q, float vr) //Calculating the Force due to collision due to wall using Linear Viscoelastic Model
{
    float F = -1*Q.k*c - 1*Q.eta*vr;
    return F;
}


int main()
{
    simDom box = simDom(0,0,0,10,10,10); // Defining the Box

    particle P = particle(3, 10000, 1000, 5, 5, 5, 1, 1, 0, 0.5, 0.25, 0.35, 0.75166); // Defining the Particle

    quaternion temp;
    sphere curr;


    ofstream file1, file2;
    file1.open("position1.txt");
    file2.open("ForceTime1.csv");

    float dt=0.00001; // One Time Step Value

    float oxy1,oxy2,oyz1,oyz2,ozx1,ozx2;
    float contact[6];
    float vr=0.0;

    for(int i=0;i<1000000;i++)
    {
        //1st Sphere
        P.A.rs = temp.product(P.q, P.A.rb);
        P.A.rs = temp.product(P.A.rs, temp.T(P.q));
        P.A.rs = temp.sum(P.rcm, P.A.rs);
        oxy1 = (P.A.rs.vz - box.z1);
        oxy2 = (-P.A.rs.vz + box.z2);
        oyz1 = (P.A.rs.vx - box.x1);
        oyz2 = (-P.A.rs.vx + box.x2);
        ozx1 = (P.A.rs.vy - box.y1);
        ozx2 = (-P.A.rs.vy + box.y2);
        contact[0] = oxy1;
        contact[1] = oxy2;
        contact[2] = oyz1;
        contact[3] = oyz2;
        contact[4] = ozx1;
        contact[5] = ozx2;
        float Farr[6];
        for(int i=0;i<6;i++)
        {
            if(i==0 && i==1)
                vr = P.vcm.vz;
            else if(i==2 && i==3)
                vr = P.vcm.vx;
            else
                vr = P.vcm.vy;

            if(contact[i]>=P.A.R)
            {
                contact[i] = 0;
            }
            else
            {
                Farr[i] = calculateForce(contact[i], P, vr);
            }
        
        
        }
        float F[3];
        F[0] = -Farr[2] + Farr[3];
        F[1] = -Farr[4] + Farr[5];
        F[2] = -Farr[0] + Farr[1];
        P.A.F.vx = F[0];
        P.A.F.vy = F[1];
        P.A.F.vz = F[2];

        //2nd Sphere
        P.B.rs = temp.product(P.q, P.B.rb);
        P.B.rs = temp.product(P.B.rs, temp.T(P.q));
        P.B.rs = temp.sum(P.rcm, P.B.rs);
        oxy1 = (P.B.rs.vz - box.z1);
        oxy2 = (-P.B.rs.vz + box.z2);
        oyz1 = (P.B.rs.vx - box.x1);
        oyz2 = (-P.B.rs.vx + box.x2);
        ozx1 = (P.B.rs.vy - box.y1);
        ozx2 = (-P.B.rs.vy + box.y2);
        contact[0] = oxy1;
        contact[1] = oxy2;
        contact[2] = oyz1;
        contact[3] = oyz2;
        contact[4] = ozx1;
        contact[5] = ozx2;        

        for(int i=0;i<6;i++)
        {
            Farr[i] = 0;

            if(i==0 && i==1)
                vr = P.vcm.vz;
            else if(i==2 && i==3)
                vr = P.vcm.vx;
            else
                vr = P.vcm.vy;

            if(contact[i]>=P.B.R)
            {
                contact[i] = 0;
            }
            else
            {
                Farr[i] = calculateForce(contact[i], P, vr);
            }

        
        }
        F[0] = -Farr[2] + Farr[3];
        F[1] = -Farr[4] + Farr[5];
        F[2] = -Farr[0] + Farr[1];
        P.B.F.vx = F[0];
        P.B.F.vy = F[1];
        P.B.F.vz = F[2];


        //3rd Sphere
        P.C.rs = temp.product(P.q, P.C.rb);
        P.C.rs = temp.product(P.C.rs, temp.T(P.q));
        P.C.rs = temp.sum(P.rcm, P.C.rs);
        oxy1 = (P.C.rs.vz - box.z1);
        oxy2 = (-P.C.rs.vz + box.z2);
        oyz1 = (P.C.rs.vx - box.x1);
        oyz2 = (-P.C.rs.vx + box.x2);
        ozx1 = (P.C.rs.vy - box.y1);
        ozx2 = (-P.C.rs.vy + box.y2);
        contact[0] = oxy1;
        contact[1] = oxy2;
        contact[2] = oyz1;
        contact[3] = oyz2;
        contact[4] = ozx1;
        contact[5] = ozx2;
        for(int i=0;i<6;i++)
        {
        Farr[i] = 0;

            if(i==0 && i==1)
                vr = P.vcm.vz;
            else if(i==2 && i==3)
                vr = P.vcm.vx;
            else
                vr = P.vcm.vy;

            if(contact[i]>=P.C.R)
            {
                contact[i] = 0;
            }
            else
            {
                Farr[i] = calculateForce(contact[i], P, vr);
            }

        }
        F[0] = -Farr[2] + Farr[3];
        F[1] = -Farr[4] + Farr[5];
        F[2] = -Farr[0] + Farr[1];
        P.C.F.vx = F[0];
        P.C.F.vy = F[1];
        P.C.F.vz = F[2];

        P.computeNetF();
        P.computeNetT();
        

        // Calculating the Linear Acceleration
        P.acm = temp.product(quaternion(1/P.m,0,0,0), P.C.F);
        P.acm.s = 0;

        // Calculating the Angular Acceleration using Euler's Equations
        P.dw.vx = (P.T.vx - (P.I.vz - P.I.vy)*P.wb.vy*P.wb.vz)/(P.I.vx);
        P.dw.vy = (P.T.vy - (P.I.vx - P.I.vz)*P.wb.vx*P.wb.vz)/(P.I.vy);
        P.dw.vz = (P.T.vz - (P.I.vy - P.I.vx)*P.wb.vx*P.wb.vy)/(P.I.vz);

        // Calculating the Rate of change of Orientation
        P.dq = temp.product(quaternion(0.5,0,0,0), temp.product(P.q, P.wb));
        // New Velocity
        P.vcm = temp.sum(P.vcm, temp.product(P.acm, quaternion(dt,0,0,0)));
        P.vcm.s = 0;
        // New Position
        P.rcm = temp.sum(P.rcm, temp.product(P.vcm, quaternion(dt,0,0,0)));
        P.rcm.s = 0;
        // New Orientation
        P.q = temp.sum(P.q, temp.product(quaternion(dt,0,0,0), P.dq));
        // New Angular Velocity
        P.wb = temp.sum(P.wb, temp.product(P.dw, quaternion(dt,0,0,0)));
        P.wb.s = 0;

        //file <<i*dt<<","<<P.rcm.vx <<","<<P.F.vx<<endl;
        if(i%5 == 0)
        {
            file1<<i*dt<<" "<<P.A.rs.vx<<" "<<P.A.rs.vy<<" "<<P.A.rs.vz<<" "<<P.B.rs.vx<<" "<<P.B.rs.vy<<" "<<P.B.rs.vz<<" "<<P.C.rs.vx<<" "<<P.C.rs.vy<<" "<<P.C.rs.vz<<" "<<P.rcm.vx<<" "<<P.rcm.vy<<" "<<P.rcm.vz<<endl;
            file2<<i*dt<<","<<P.rcm.vx<<","<<P.F.vx<<","<<P.vcm.vx<<endl;
        }


    }
    file1.close();
    file2.close();
    return 0;
}

