#include <thrust/device_vector.h>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#include <fstream>
#include <stdio.h>
//#include <string.h>
#include <sstream>

typedef thrust::device_vector< int> int_vec_type;
typedef double value_type;
struct er3bp_system
{
    struct er3bp_functor
    {
        double m_t;
        double m_mu;
        double m_ecc;
        er3bp_functor( double t, double mu, double ecc ) : m_t(t),m_mu(mu),m_ecc(ecc){}

        template< class T >
        __host__ __device__
        void operator()( T tpl ) const
        {
            value_type c=thrust::get< 8 >(tpl);
            value_type x = thrust::get< 0 >( tpl );
            value_type y = thrust::get< 1 >( tpl );
            value_type vx = thrust::get< 2 >( tpl );
            value_type vy = thrust::get< 3 >( tpl );
            thrust::get< 4 >(tpl)=vx; // set tuple's fifth element to vx
            thrust::get< 5 >(tpl)=vy;
            thrust::get< 6 >(tpl)=2*vy+(x-((1-m_mu)*(x+m_mu))/pow((x+m_mu)*(x+m_mu)+y*y,1.5)-
                                        (m_mu*(x-1+m_mu))/pow((x-1+m_mu)*(x-1+m_mu)+y*y,1.5))/(1+m_ecc*cos(m_t));
            thrust::get< 7 >(tpl)=-2*vx+(y-(1-m_mu)*y/pow((x+m_mu)*(x+m_mu)+y*y,1.5)-
                                         m_mu*y/pow((x-1+m_mu)*(x-1+m_mu)+y*y,1.5))/(1+m_ecc*cos(m_t));
            //printf("%.12f\n",m_t);
            if (c==0){
            if ( (sqrt( (x-m_mu)*(x-m_mu)+y*y)<0.015) || (sqrt( (x-1+m_mu)*(x-1+m_mu)+y*y)<0.0044)) //FIXME aggiustare formula e mettere numero adimensionale
            {
                thrust::get< 8 >(tpl)=1;
                //printf("collision!");
            }
            }
            }
    };


    er3bp_system( size_t N, int_vec_type &coll ): m_N( N ), m_coll(coll) { }

    template< class State , class Deriv >
    void operator()( State &x , Deriv &dxdt , value_type t ) const
    {
        thrust::for_each(
        // Create each tuple
                thrust::make_zip_iterator( thrust::make_tuple(
                        boost::begin( x ) ,             // x
                        boost::begin( x ) + m_N ,       // y
                        boost::begin( x ) + 2 * m_N ,   // vx
                        boost::begin( x ) + 3 * m_N,    // vy
                        boost::begin( dxdt ) ,
                        boost::begin( dxdt ) + m_N ,
                        boost::begin( dxdt ) + 2 * m_N,
                        boost::begin( dxdt ) + 3 * m_N
                        ,boost::begin(m_coll) ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple(
                        boost::begin( x ) + m_N ,
                        boost::begin( x ) + 2 * m_N ,
                        boost::begin( x ) + 3 * m_N ,
                        boost::begin( x ) + 4 * m_N ,
                        boost::begin( dxdt ) + m_N ,
                        boost::begin( dxdt ) + 2 * m_N ,
                        boost::begin( dxdt ) + 3 * m_N ,
                        boost::begin( dxdt ) + 4 * m_N
                        ,boost::begin(m_coll)+m_N ) ) ,
                er3bp_functor(t, mu, ecc) );
                //printf("%.12f\n",t);
    }

size_t m_N;
int_vec_type &m_coll;
};
