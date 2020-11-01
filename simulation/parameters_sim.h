#ifndef P_SIM_H
#define P_SIM_H

#include <QObject>
#include <QDebug>

#include <Eigen/Core>
#include <iostream>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/archive/basic_binary_iarchive.hpp>
#include <boost/archive/basic_binary_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

// variables related to the formulation of the model

namespace icy { class SimParams; }

class icy::SimParams : public QObject
{
    Q_OBJECT

    // general
    //Q_PROPERTY(bool s_SaveResult MEMBER SaveResult NOTIFY propertyChanged)
    Q_PROPERTY(int s_MaxSteps MEMBER MaxSteps NOTIFY propertyChanged)

    // integration
    Q_PROPERTY(double in_InitialTimeStep MEMBER InitialTimeStep NOTIFY propertyChanged)
    Q_PROPERTY(double in_HHTalpha READ getHHTalpha WRITE setHHTalpha)
    Q_PROPERTY(double in_NewmarkBeta READ getNewmarkBeta)
    Q_PROPERTY(double in_NewmarkGamma READ getNewmarkGamma)

    Q_PROPERTY(double in_Damping MEMBER Damping NOTIFY propertyChanged)

    Q_PROPERTY(double in_ConvergenceEpsilon MEMBER ConvergenceEpsilon NOTIFY propertyChanged)
    Q_PROPERTY(double in_ConvergenceCutoff MEMBER ConvergenceCutoff NOTIFY propertyChanged)
    Q_PROPERTY(int in_IterationsMax MEMBER IterationsMax NOTIFY propertyChanged)
    Q_PROPERTY(int in_IterationsMin MEMBER IterationsMin NOTIFY propertyChanged)

    // material parameters and physical constants
    Q_PROPERTY(double p_Gravity_z MEMBER gravity NOTIFY propertyChanged)
    Q_PROPERTY(double p_WaterDensity MEMBER WaterDensity NOTIFY propertyChanged)
    Q_PROPERTY(double p_IceDensity MEMBER IceDensity NOTIFY propertyChanged)
    Q_PROPERTY(double p_YoungsModulus MEMBER YoungsModulus NOTIFY propertyChanged)
    Q_PROPERTY(double p_PoissonsRatio MEMBER PoissonsRatio NOTIFY propertyChanged)
    Q_PROPERTY(double p_Thickness MEMBER Thickness NOTIFY propertyChanged)

    // fracture
    Q_PROPERTY(double f_normal_traction_threshold MEMBER normal_traction_threshold NOTIFY propertyChanged)
    Q_PROPERTY(double f_epsilon MEMBER fracture_epsilon NOTIFY propertyChanged)
    Q_PROPERTY(double f_CharacteristicLengthMax MEMBER CharacteristicLengthMax NOTIFY propertyChanged)
    Q_PROPERTY(int f_substep_radius MEMBER substep_radius NOTIFY propertyChanged)
    Q_PROPERTY(int f_substep_radius2 MEMBER substep_radius2 NOTIFY propertyChanged)
    Q_PROPERTY(bool f_enable MEMBER fracture_enable NOTIFY propertyChanged)
    Q_PROPERTY(int f_max_substeps MEMBER fracture_max_substeps NOTIFY propertyChanged)
    Q_PROPERTY(int f_substep_iterations MEMBER substep_iterations NOTIFY propertyChanged)
    Q_PROPERTY(double f_weakening MEMBER weakening_coeff NOTIFY propertyChanged)
    Q_PROPERTY(double f_temporal_attenuation MEMBER temporal_attenuation NOTIFY propertyChanged)
    Q_PROPERTY(double f_wave_height MEMBER wave_height NOTIFY propertyChanged)
    Q_PROPERTY(double f_wave_start_location MEMBER wave_start_location NOTIFY propertyChanged)
    Q_PROPERTY(double f_substepping_timestep_factor MEMBER substepping_timestep_factor NOTIFY propertyChanged)

    // load type
    Q_PROPERTY(int t_Load READ getLoadType)

public:
    double InitialTimeStep, HHTalpha, NewmarkBeta, NewmarkGamma;
    double ConvergenceEpsilon,  ConvergenceCutoff;
    int MaxSteps, IterationsMax, IterationsMin;
    double gravity, WaterDensity, IceDensity, PoissonsRatio, YoungsModulus;
    double Damping, Thickness;
    bool SaveResult;

    int loadType;       // for testing

    // fracture
    float normal_traction_threshold;
    float fracture_epsilon;
    double CharacteristicLengthMax; // initial meshing
    int substep_radius, substep_radius2; // how many neighbor levels involved in local substep
    bool fracture_enable;
    int fracture_max_substeps;
    int substep_iterations;
    float weakening_coeff;
    double temporal_attenuation;
    double wave_height;
    double wave_start_location;
    double substepping_timestep_factor;
    double cutoff_coefficient; // used to reduce the computation cost of max_traction

    // other

    SimParams() { Reset(); }

    void Reset() {
        SaveResult = false;
        MaxSteps = 2000;

        // integration
        InitialTimeStep = 0.02;

        ConvergenceEpsilon = 0.01;
        ConvergenceCutoff = 1E-8;
        IterationsMax = 10;
        IterationsMin = 2;

        // material parameters and physical constants
        gravity = 9.81;
        WaterDensity = 997;
        IceDensity = 910;

        Damping = 0.1;
        Thickness = 0.1;

        loadType = 0;

        PoissonsRatio = 0.3;
        YoungsModulus = 3.7e9;

        setHHTalpha(0.3);

        normal_traction_threshold = 1e5;
        fracture_epsilon = 0.1;
        CharacteristicLengthMax = 0.5;
        substep_radius = 10;
        substep_radius2 = 150;
        fracture_enable = true;
        fracture_max_substeps=1000;
        substep_iterations = 2;
        weakening_coeff = 0.75;
        temporal_attenuation = 0.2;
        wave_height = 0.1;
        substepping_timestep_factor = 0.001;
        wave_start_location = 10;
        cutoff_coefficient = 0.4;

        emit propertyChanged();
    }

    // integration
    double getHHTalpha() {return HHTalpha;}
    void setHHTalpha(double alpha) {
        HHTalpha = alpha;
        NewmarkBeta = (1+alpha)*(1+alpha)*0.25;
        NewmarkGamma = 0.5+alpha;
        emit propertyChanged();
    }
    double getNewmarkBeta() {return NewmarkBeta;}
    double getNewmarkGamma() {return NewmarkGamma;}
    int getLoadType() {return loadType;}

    // serialization
    static unsigned const buffer_size = 400;
    char serialization_buffer[buffer_size];

    void Serialize()
    {
        boost::iostreams::basic_array_sink<char> sr(serialization_buffer, buffer_size);
        boost::iostreams::stream< boost::iostreams::basic_array_sink<char> > source(sr);
        boost::archive::binary_oarchive oa(source);

        oa << SaveResult;
        oa << MaxSteps;
        oa << InitialTimeStep;
        oa << HHTalpha;
        oa << ConvergenceEpsilon;
        oa << ConvergenceCutoff;
        oa << IterationsMax;
        oa << IterationsMin;
        oa << gravity;
        oa << WaterDensity;
        oa << IceDensity;
        oa << PoissonsRatio;
        oa << Damping;
        oa << Thickness;
        oa << loadType;
        oa << YoungsModulus;
        oa << normal_traction_threshold;
        oa << fracture_epsilon;
        oa << CharacteristicLengthMax;
        oa << substep_radius;
        oa << fracture_enable;
        oa << fracture_max_substeps;
        oa << substep_iterations;
        oa << weakening_coeff;
        oa << temporal_attenuation;
        oa << wave_height;
        oa << wave_start_location;
        oa << substepping_timestep_factor;
        oa << substep_radius2;
    }

    void Deserialize()
    {
        boost::iostreams::basic_array_source<char> src(serialization_buffer, buffer_size);
        boost::iostreams::stream< boost::iostreams::basic_array_source<char> > source(src);
        boost::archive::binary_iarchive oa(source);
        oa >> SaveResult;
        oa >> MaxSteps;
        oa >> InitialTimeStep;
        oa >> HHTalpha;
        oa >> ConvergenceEpsilon;
        oa >> ConvergenceCutoff;
        oa >> IterationsMax;
        oa >> IterationsMin;
        oa >> gravity;
        oa >> WaterDensity;
        oa >> IceDensity;
        oa >> PoissonsRatio;
        oa >> Damping;
        oa >> Thickness;
        oa >> loadType;
        oa >> YoungsModulus;
        oa >> normal_traction_threshold;
        oa >> fracture_epsilon;
        oa >> CharacteristicLengthMax;
        oa >> substep_radius;
        oa >> fracture_enable;
        oa >> fracture_max_substeps;
        oa >> substep_iterations;
        oa >> weakening_coeff;
        oa >> temporal_attenuation;
        oa >> wave_height;
        oa >> wave_start_location;
        oa >> substepping_timestep_factor;
        oa >> substep_radius2;
    }


signals:
    void propertyChanged();
};


#endif
