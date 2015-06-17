#pragma once

#include "base.h"

#include <vector>
#include <memory>
#include <algorithm>

class CompositePotential : public Potential {

public:
    std::vector< std::unique_ptr< Potential > > potentials;

    CompositePotential():potentials()
    {}

    std::string get_name() const
    {
        return "CompositePotential";
    }

    inline size_t num_potentials() const
    {
        return potentials.size();
    }

    std::unique_ptr< Potential > convert_to_single_potential()
    {
        return std::move(potentials[0]);
    }

    double operator()(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return (*p)(x,y,z,t);} );
    }

    double ddx(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->ddx(x,y,z,t);} );
    }

    double ddy(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->ddy(x,y,z,t);} );
    }

    double ddz(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->ddz(x,y,z,t);} );
    }

    double d2dr2(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->d2dr2(x,y,z,t);} );
    }

    double d2dx2(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->d2dx2(x,y,z,t);} );
    }

    double d2dxdy(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->d2dxdy(x,y,z,t);} );
    }

    double d2dxdz(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->d2dxdz(x,y,z,t);} );
    }

    double d2dy2(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->d2dy2(x,y,z,t);} );
    }

    double d2dydz(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->d2dydz(x,y,z,t);} );
    }

    double d2dz2(double x, double y, double z, double t) const {
        return xform_accumulate( 0.0, [=](const std::unique_ptr<Potential>& p)
                                                {return p->d2dz2(x,y,z,t);} );
    }

    void acc(double x, double y, double z, double t, Vec3D& acc) const {
        acc.x = 0;
        acc.y = 0;
        acc.z = 0;

        Vec3D toadd;
        for( const auto& p : potentials ) {
            p->acc(x, y, z, t, toadd);
            acc += toadd;
        }
    }

private:
    template <typename T, typename K>
    T xform_accumulate(T init, const K xformer) const
    {
        return std::accumulate( potentials.cbegin(), potentials.cend(), init,
                                [xformer](T a, const std::unique_ptr<Potential>& p)
                                         {return a + xformer(p);}
                              );
    }

};



