#ifndef FLOEVISUALIZATION_H
#define FLOEVISUALIZATION_H

#include <vector>

#include <QObject>
#include <QString>
#include <QDebug>
#include <QMetaEnum>

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>

#include <vtkLabeledDataMapper.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkIdFilter.h>
#include <vtkCellCenters.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>

#include <vtkStructuredGrid.h>

#include "geometry.h"

namespace icy { class FloeVisualization; class Node; class Element; class Edge;}

class icy::FloeVisualization : public QObject
{
    Q_OBJECT
public:
    FloeVisualization();

    // visualization options
    enum VisOpt { none, boundary, vert_force,
                  fracture_support, max_normal_traction,
                  deflection, AbsMx, Mx, My, Mxy, Mx_e, My_e, Mxy_e, Tx, Ty, Txy, Qx, Qy,
                  stx, sty, stxy, st1, st2, sbx, sby, sbxy, sb1, sb2};
    Q_ENUM(VisOpt)

    // presenting via VTK
    bool use_tentative_coordinates = true;
    bool update_minmax = true;
    bool update_arrows = true;

    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkPoints> points;

    // 2D mesh
    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkCellArray> cellArray;
    vtkNew<vtkDataSetMapper> dataSetMapper;
    vtkNew<vtkActor> actor_mesh;

    // 3D mesh
    vtkNew<vtkPoints> points_mesh3d;
    vtkNew<vtkUnstructuredGrid> ugrid_mesh3d;
    vtkNew<vtkCellArray> cellArray_mesh3d;
    vtkNew<vtkDataSetMapper> mapper_mesh3d;
    vtkNew<vtkActor> actor_mesh_mesh3d;

    // labels
    vtkNew<vtkUnstructuredGrid> ugrid_selection;
    vtkNew<vtkIntArray> edgeNumbers;
    vtkNew<vtkUnstructuredGridGeometryFilter> geometryFilter;
    vtkNew<vtkLabeledDataMapper> labledDataMapper;
    vtkNew<vtkActor2D> actor_labels;
    vtkNew<vtkIdFilter> idfilter;
    vtkNew<vtkCellCenters> cellCenters;

    // boundary
    vtkNew<vtkCellArray> cellArray_boundary;
    vtkNew<vtkUnstructuredGrid> ugrid_boundary;
    vtkNew<vtkDataSetMapper> dataSetMapper_boundary;
    vtkNew<vtkActor> actor_boundary;

    // arrows
    vtkNew<vtkCellArray> cellArray_vertices;
    vtkNew<vtkUnstructuredGrid> ugrid_vertices;
    vtkNew<vtkDoubleArray> arrowCoords;
    vtkNew<vtkArrowSource> arrowSource;
    vtkNew<vtkGlyph3D> glyph3D;
    vtkNew<vtkPolyDataMapper> mapper_arrows;
    vtkNew<vtkActor> actor_arrows;

    // surface
    static constexpr unsigned gridSizeX = 100;
    static constexpr unsigned gridSizeY = 50;
    vtkNew<vtkStructuredGrid> grid_water;
    vtkNew<vtkPoints> points_water;
    vtkNew<vtkDataSetMapper> mapper_water;
    vtkNew<vtkActor> actor_water;

    // visualizing variables
    vtkNew<vtkDoubleArray> visualized_values;
    vtkIdType selectedPointId = -1;
    // visual representation of the floe with thickness
    void UnsafeUpdateTopology(std::vector<icy::Node*> *nodes, std::vector<icy::Element*> *elems,
                              std::vector<icy::Edge*> *edges);
    void UnsafeUpdateDisplacements(std::vector<icy::Node*> *nodes, std::vector<icy::Element*> *elems);
    void UnsafeUpdateValues(std::vector<icy::Node*> *nodes,
                            std::vector<icy::Element*> *elems,
                            int option = -1);
    void UnsafeUpdateArrows(std::vector<icy::Node*> *nodes);
    void UnsafeUpdateSelection(std::vector<icy::Node*> *nodes,
                               vtkIdType selectedPoint = -1);
    void UnsafeUpdateWaterLine(double simulationTime, SimParams &prms);



    // visualization
    VisOpt VisualizingVariable = VisOpt::none;
    vtkNew<vtkLookupTable> hueLut;

    void InitializeLUT(int table);

    vtkNew<vtkLookupTable> lutAsymmetry;

    static constexpr float lutArrayTerrain[51][3] =
    {{0.54938, 0.772213, 0.848103},
     {0.54349, 0.751689, 0.811536},
     {0.5376, 0.731165, 0.77497},
     {0.53171, 0.71064, 0.738403},
     {0.525821, 0.690116, 0.701836},
     {0.523788, 0.674231, 0.670908},
     {0.526385, 0.663912, 0.646744},
     {0.528982, 0.653593, 0.622581},
     {0.531579, 0.643275, 0.598418},
     {0.534176, 0.632956, 0.574255},
     {0.542585, 0.630953, 0.560901},
     {0.551575, 0.629781, 0.548628},
     {0.560566, 0.628609, 0.536356},
     {0.569556, 0.627438, 0.524083},
     {0.579925, 0.628775, 0.515402},
     {0.592706, 0.634504, 0.513005},
     {0.605487, 0.640233, 0.510609},
     {0.618268, 0.645962, 0.508213},
     {0.631049, 0.651691, 0.505817},
     {0.644928, 0.660899, 0.509458},
     {0.65905, 0.67088, 0.514441},
     {0.673172, 0.680861, 0.519424},
     {0.687294, 0.690842, 0.524407},
     {0.701252, 0.701293, 0.530766},
     {0.71477, 0.712998, 0.540793},
     {0.728289, 0.724704, 0.55082},
     {0.741807, 0.736409, 0.560848},
     {0.755325, 0.748114, 0.570875},
     {0.767446, 0.759552, 0.583215},
     {0.779043, 0.77089, 0.596421},
     {0.79064, 0.782228, 0.609627},
     {0.802237, 0.793566, 0.622834},
     {0.813354, 0.804566, 0.636368},
     {0.822309, 0.81404, 0.651377},
     {0.831264, 0.823514, 0.666387},
     {0.84022, 0.832989, 0.681397},
     {0.849175, 0.842463, 0.696406},
     {0.8563, 0.850209, 0.711857},
     {0.862381, 0.856968, 0.727562},
     {0.868461, 0.863726, 0.743266},
     {0.874541, 0.870485, 0.75897},
     {0.880373, 0.876977, 0.774624},
     {0.883714, 0.880806, 0.789783},
     {0.887055, 0.884635, 0.804943},
     {0.890396, 0.888464, 0.820102},
     {0.893737, 0.892293, 0.835261},
     {0.895825, 0.894749, 0.849094},
     {0.896869, 0.896062, 0.86182},
     {0.897913, 0.897375, 0.874547},
     {0.898956, 0.898687, 0.887273},
     {0.9, 0.9, 0.9}};

    static constexpr float lutArrayThermometer[51][3] = {
        {0.163302, 0.119982, 0.79353},
        {0.183304, 0.162484, 0.815377},
        {0.203306, 0.204986, 0.837223},
        {0.223309, 0.247488, 0.85907},
        {0.243311, 0.28999, 0.880917},
        {0.269511, 0.336207, 0.897377},
        {0.303148, 0.386882, 0.907374},
        {0.336786, 0.437557, 0.917372},
        {0.370423, 0.488231, 0.927369},
        {0.404061, 0.538906, 0.937366},
        {0.440238, 0.581486, 0.94075},
        {0.476669, 0.623257, 0.943472},
        {0.5131, 0.665028, 0.946195},
        {0.549532, 0.706799, 0.948918},
        {0.584528, 0.743128, 0.948889},
        {0.617013, 0.769936, 0.944046},
        {0.649498, 0.796744, 0.939202},
        {0.681983, 0.823552, 0.934359},
        {0.714468, 0.85036, 0.929516},
        {0.74029, 0.863925, 0.917779},
        {0.764631, 0.874548, 0.904511},
        {0.788973, 0.88517, 0.891243},
        {0.813314, 0.895793, 0.877975},
        {0.834823, 0.902115, 0.862594},
        {0.848779, 0.896972, 0.841576},
        {0.862735, 0.891829, 0.820558},
        {0.87669, 0.886685, 0.79954},
        {0.890646, 0.881542, 0.778522},
        {0.89662, 0.865623, 0.753608},
        {0.8996, 0.845665, 0.727233},
        {0.90258, 0.825706, 0.700858},
        {0.905561, 0.805747, 0.674483},
        {0.906659, 0.78341, 0.647684},
        {0.899291, 0.750375, 0.618978},
        {0.891922, 0.717339, 0.590273},
        {0.884553, 0.684303, 0.561568},
        {0.877185, 0.651267, 0.532862},
        {0.863955, 0.612218, 0.504161},
        {0.847377, 0.569734, 0.475461},
        {0.830799, 0.527249, 0.446761},
        {0.814221, 0.484764, 0.418062},
        {0.796881, 0.442068, 0.389497},
        {0.771921, 0.397256, 0.362281},
        {0.746961, 0.352443, 0.335064},
        {0.722001, 0.307631, 0.307848},
        {0.697041, 0.262818, 0.280632},
        {0.667501, 0.223593, 0.256072},
        {0.634146, 0.189023, 0.233727},
        {0.600791, 0.154453, 0.211381},
        {0.567436, 0.119883, 0.189036},
        {0.534081, 0.0853132, 0.16669}};


    static constexpr float lutArrayTemperature[51][3] =
{{0.165698, 0.282261, 0.936187}, {0.200314, 0.334793,
  0.938464}, {0.234929, 0.387325, 0.940741}, {0.269545, 0.439857,
  0.943017}, {0.30416, 0.492389, 0.945294}, {0.338776, 0.544921,
  0.947571}, {0.373072, 0.590873, 0.950256}, {0.406968, 0.6286,
  0.95345}, {0.440864, 0.666328, 0.956645}, {0.47476, 0.704056,
  0.959839}, {0.508657, 0.741783, 0.963034}, {0.542553, 0.779511,
  0.966228}, {0.573124, 0.806257, 0.969024}, {0.603279, 0.83163,
  0.971771}, {0.633434, 0.857004, 0.974518}, {0.663589, 0.882377,
  0.977264}, {0.693745, 0.90775, 0.980011}, {0.722982, 0.928921,
  0.982462}, {0.750385, 0.941688, 0.984321}, {0.777788, 0.954454,
  0.986181}, {0.805191, 0.96722, 0.98804}, {0.832594, 0.979987,
  0.989899}, {0.859997, 0.992753, 0.991759}, {0.884446, 0.995483,
  0.968897}, {0.908051, 0.995346, 0.938973}, {0.931655, 0.995208,
  0.909049}, {0.955259, 0.995071, 0.879124}, {0.978864, 0.994934,
  0.8492}, {0.996605, 0.992863, 0.819085}, {0.993822, 0.984024,
  0.788303}, {0.99104, 0.975186, 0.757522}, {0.988257, 0.966347,
  0.72674}, {0.985475, 0.957509, 0.695959}, {0.982692, 0.94867,
  0.665177}, {0.977149, 0.930487, 0.633887}, {0.970225, 0.907633,
  0.602342}, {0.9633, 0.884778, 0.570797}, {0.956376, 0.861923,
  0.539252}, {0.949452, 0.839068, 0.507707}, {0.94233, 0.814829,
  0.476285}, {0.933626, 0.779513, 0.445847}, {0.924922, 0.744197,
  0.415409}, {0.916218, 0.708881, 0.384971}, {0.907514, 0.673565,
  0.354533}, {0.89881, 0.638249, 0.324095}, {0.889975, 0.597251,
  0.298625}, {0.881035, 0.551706, 0.27713}, {0.872095, 0.506162,
  0.255635}, {0.863154, 0.460617, 0.23414}, {0.854214, 0.415073,
  0.212645}, {0.845274, 0.369528, 0.19115}};

    static constexpr float lutArrayTemperatureAdj[51][3] =
{{0.770938, 0.951263, 0.985716}, {0.788065, 0.959241,
  0.986878}, {0.805191, 0.96722, 0.98804}, {0.822318, 0.975199,
  0.989202}, {0.839445, 0.983178, 0.990364}, {0.856572, 0.991157,
  0.991526}, {0.872644, 0.995552, 0.98386}, {0.887397, 0.995466,
  0.965157}, {0.902149, 0.99538, 0.946454}, {0.916902, 0.995294,
  0.927751}, {0.931655, 0.995208, 0.909049}, {0.946408, 0.995123,
  0.890346}, {0.961161, 0.995037, 0.871643}, {0.975913, 0.994951,
  0.85294}, {0.990666, 0.994865, 0.834237}, {0.996257, 0.991758,
  0.815237}, {0.994518, 0.986234, 0.795999}, {0.992779, 0.98071,
  0.77676}, {0.99104, 0.975186, 0.757522}, {0.989301, 0.969662,
  0.738283}, {0.987562, 0.964138, 0.719045}, {0.985823, 0.958614,
  0.699807}, {0.984084, 0.953089, 0.680568}, {0.982345, 0.947565,
  0.66133}, {0.97888, 0.936201, 0.641773}, {0.974552, 0.921917,
  0.622058}, {0.970225, 0.907633, 0.602342}, {0.965897, 0.893348,
  0.582626}, {0.961569, 0.879064, 0.562911}, {0.957242, 0.86478,
  0.543195}, {0.952914, 0.850496, 0.52348}, {0.948586, 0.836212,
  0.503764}, {0.944259, 0.821927, 0.484048}, {0.939066, 0.801586,
  0.464871}, {0.933626, 0.779513, 0.445847}, {0.928186, 0.757441,
  0.426823}, {0.922746, 0.735368, 0.4078}, {0.917306, 0.713296,
  0.388776}, {0.911866, 0.691223, 0.369752}, {0.906426, 0.669151,
  0.350728}, {0.900986, 0.647078, 0.331704}, {0.895546, 0.625006,
  0.312681}, {0.889975, 0.597251, 0.298625}, {0.884388, 0.568785,
  0.285191}, {0.8788, 0.54032, 0.271756}, {0.873212, 0.511855,
  0.258322}, {0.867625, 0.483389, 0.244888}, {0.862037, 0.454924,
  0.231453}, {0.856449, 0.426459, 0.218019}, {0.850862, 0.397993,
  0.204584}, {0.845274, 0.369528, 0.19115}};

    static constexpr float lutArrayBands[6][3] =
{{0.684154, 0.875059, 0.95523}, {0.882734, 0.823682,
                                 0.375036},{0.601724, 0.748252,
  0.492262}, {0.745543, 0.423584, 0.52455}, {0.756931, 0.53206,
  0.359658}, {0.384161, 0.395036, 0.599407} };

};

#endif // FLOEVISUALIZATION_H
