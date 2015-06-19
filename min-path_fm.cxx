////program for itkSpeedFunctionToPathFilter (Fast Marching Minimal Path Extraction)
//01: based on template.cxx and vo-img2v-skel_01.cxx


#include <string>
#include <sstream>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkSpeedFunctionToPathFilter.h>
#include <itkGradientDescentOptimizer.h>
#include <itkPathIterator.h>
#include <itkPolyLineParametricPath.h>
#include <itkImageFileWriter.h>

#include <itkMesh.h>
#include <itkLineCell.h>
#include <itkMeshFileWriter.h>


template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    const char offset= 5;
    if( argc != offset + 3*Dimension){
        fprintf(stderr, "%d + 3*Dimension = %d parameters are needed!\n", offset-1, offset-1 + 3*Dimension);
        return EXIT_FAILURE;
        }

    char* outPrefix= argv[2];
    std::stringstream sss;

    typedef double   SpeedPixelType;
    typedef uint8_t  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<SpeedPixelType, Dimension>  SpeedImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    reader->ReleaseDataFlagOn();
    FilterWatcher watcherI(reader);
    watcherI.QuietOn();
    watcherI.ReportTimeOn();
    try{
        reader->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    // scale image values to be in [0; 1]
    typedef typename itk::RescaleIntensityImageFilter<InputImageType, SpeedImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(reader->GetOutput());
    rescaleFilter->SetOutputMinimum(0.0);
    rescaleFilter->SetOutputMaximum(1.0);
    FilterWatcher watcherR(rescaleFilter);
    try {
        rescaleFilter->Update();
        }
    catch (itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    const typename SpeedImageType::Pointer& input= rescaleFilter->GetOutput();

    typedef itk::PolyLineParametricPath<Dimension> PathType;
    typedef itk::SpeedFunctionToPathFilter<SpeedImageType, PathType> PathFilterType;
    typedef typename PathFilterType::CostFunctionType::CoordRepType CoordRepType;

    // Create interpolator
    typedef itk::LinearInterpolateImageFunction<SpeedImageType, CoordRepType> InterpolatorType;
    typename InterpolatorType::Pointer interp = InterpolatorType::New();

    // Create cost function
    typename PathFilterType::CostFunctionType::Pointer cost = PathFilterType::CostFunctionType::New();
    cost->SetInterpolator(interp);

    // Create optimizer
    typedef itk::GradientDescentOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetNumberOfIterations(atoi(argv[4]));

    // Create path filter
    typename PathFilterType::Pointer pathFilter = PathFilterType::New();
    pathFilter->SetInput(input); //needs the image values to be scaled to [0; 1]
    pathFilter->SetCostFunction(cost);
    pathFilter->SetOptimizer(optimizer);
    pathFilter->SetTerminationValue(2.0);

    // Setup path points
    typename SpeedImageType::IndexType start, end, way;

    for(int j= 0; j < Dimension; j++){
        start[j]= atoi(argv[j+offset]) - 1;
        }

    for(int j= 0; j < Dimension; j++){
        end[j]= atoi(argv[j+Dimension+offset]) - 1;
        }

    for(int j= 0; j < Dimension; j++){
        way[j]= atoi(argv[j+2*Dimension+offset]) - 1;
        }

    typename PathFilterType::PointType startP, endP, wayP;

    input->TransformIndexToPhysicalPoint(start, startP);
    input->TransformIndexToPhysicalPoint(end, endP);
    input->TransformIndexToPhysicalPoint(way, wayP);

    // Add path information
    typename PathFilterType::PathInfo info;
    info.SetStartPoint(startP);
    info.SetEndPoint(endP);
    info.AddWayPoint(wayP);
    pathFilter->AddPathInfo(info);


    // Compute the path
    FilterWatcher watcher(pathFilter); //filter reports no progress so far
    try {
        pathFilter->Update();
        }
    catch (itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    // Allocate output image
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

    typename OutputImageType::Pointer output = OutputImageType::New();
    output->SetRegions(input->GetLargestPossibleRegion());
    output->SetSpacing(input->GetSpacing());
    output->SetOrigin(input->GetOrigin());
    output->Allocate();
    output->FillBuffer(itk::NumericTraits<OutputPixelType>::Zero);

    // Rasterize path
    typedef itk::PathIterator<OutputImageType, PathType> PathIteratorType;
    std::cout << "# of paths: " << pathFilter->GetNumberOfOutputs() << std::endl;
    for (unsigned int i= 0; i < pathFilter->GetNumberOfOutputs(); i++){

        // Get the path
        typename PathType::Pointer path = pathFilter->GetOutput(i);

        // Check path is valid
        if (path->GetVertexList()->Size() == 0){
            std::cerr << "WARNING: Path " << (i+1) << " contains no points!" << std::endl;
            continue;
            }

        printf("Path %3d contains %6d points.\n", i, path->GetVertexList()->Size());

        // Iterate path and convert to image
        PathIteratorType it(output, path);
        for (it.GoToBegin(); !it.IsAtEnd(); ++it){
            it.Set(itk::NumericTraits<OutputPixelType>::max()/2);
            }
        }

    output->SetPixel(start,itk::NumericTraits<OutputPixelType>::max());
    output->SetPixel(end,itk::NumericTraits<OutputPixelType>::max());
    output->SetPixel(way,itk::NumericTraits<OutputPixelType>::max());


    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    sss << outPrefix << ".mha";
    writer->SetFileName(sss.str().c_str());
    writer->SetInput(output);
    writer->SetUseCompression(atoi(argv[3]));
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    //// create mesh to save in a VTK-file
    typedef typename itk::Mesh<float, Dimension>  MeshType;
    typename MeshType::Pointer  mesh = MeshType::New();

    for (unsigned int i=0; i < pathFilter->GetNumberOfOutputs(); i++){

        // Get the path
        typename PathType::Pointer path = pathFilter->GetOutput(i);
        const typename PathType::VertexListType *vertexList = path->GetVertexList();

        for(unsigned int k = 0; k < vertexList->Size(); k++){
            mesh->SetPoint(k, vertexList->GetElement(k));
            }
        }

    std::cout << "# of mesh points: " << mesh->GetNumberOfPoints() << std::endl;

    //// create connecting lines
    typedef typename MeshType::CellType CellType;
    typedef typename itk::LineCell<CellType> LineType;
    typedef typename CellType::CellAutoPointer  CellAutoPointer;

    //// from: http://www.itk.org/Doxygen/html/Examples_2DataRepresentation_2Mesh_2Mesh3_8cxx-example.html
    const unsigned int numberOfCells = mesh->GetNumberOfPoints() - 1;
    typename CellType::CellAutoPointer line;
    for(size_t cellId=0; cellId < numberOfCells; cellId++){
        line.TakeOwnership(new LineType);
        line->SetPointId(0, cellId);
        line->SetPointId(1, cellId+1);
        mesh->SetCell(cellId, line);
        }

    std::cout << "# of mesh cells: " << mesh->GetNumberOfCells() << std::endl;

    typedef typename itk::MeshFileWriter<MeshType> MeshWriterType;
    typename MeshWriterType::Pointer mwriter = MeshWriterType::New();

    FilterWatcher watcherMO(mwriter);
    sss.str(""); sss << outPrefix << ".vtk"; //vtp not supported as of itk-4.8
    mwriter->SetFileName(sss.str().c_str());
    mwriter->SetInput(mesh);
    mwriter->Update();


    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= 0;
    switch (dimensionType){
    case 2:
        res= DoIt<InputComponentType, InputPixelType, 2>(argc, argv);
        break;
    case 3:
        res= DoIt<InputComponentType, InputPixelType, 3>(argc, argv);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= 0;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

int dispatch_cT(itk::ImageIOBase::IOComponentType componentType, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= 0;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::FLOAT:{        // float32
        typedef float InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::DOUBLE:{       // float64
        typedef double InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    default:
        std::cerr << "unknown component type" << std::endl;
        break;
        }//switch
    return res;
    }


////from http://itk-users.7.n7.nabble.com/Pad-image-with-0-but-keep-its-type-what-ever-it-is-td27442.html
//namespace itk{
  // Description:
  // Get the PixelType and ComponentType from fileName

void GetImageType (std::string fileName,
    itk::ImageIOBase::IOPixelType &pixelType,
    itk::ImageIOBase::IOComponentType &componentType,
    size_t &dimensionType
    ){
    typedef itk::Image<char, 1> ImageType; //template initialization parameters need to be given but can be arbitrary here
    itk::ImageFileReader<ImageType>::Pointer imageReader= itk::ImageFileReader<ImageType>::New();
    imageReader->SetFileName(fileName.c_str());
    imageReader->UpdateOutputInformation();

    pixelType = imageReader->GetImageIO()->GetPixelType();
    componentType = imageReader->GetImageIO()->GetComponentType();
    dimensionType= imageReader->GetImageIO()->GetNumberOfDimensions();

    std::cerr << std::endl << "dimensions: " << dimensionType << std::endl;
    std::cerr << "component type: " << imageReader->GetImageIO()->GetComponentTypeAsString(componentType) << std::endl;
    std::cerr << "component size: " << imageReader->GetImageIO()->GetComponentSize() << std::endl;
    std::cerr << "pixel type (string): " << imageReader->GetImageIO()->GetPixelTypeAsString(imageReader->GetImageIO()->GetPixelType()) << std::endl;
    std::cerr << "pixel type: " << pixelType << std::endl << std::endl;

    }



int main(int argc, char *argv[]){
    if ( argc < 5 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image_Base"
                  << " compress"
                  << " iterations"
                  << " start-point..."
                  << " end-point..."
                  << " way-point..."
                  << std::endl;
        std::cerr << " Point coordinates are expected in voxel units (starting with 1)!" << std::endl;

        return EXIT_FAILURE;
        }

    itk::ImageIOBase::IOPixelType pixelType;
    typename itk::ImageIOBase::IOComponentType componentType;
    size_t dimensionType;


    try {
        GetImageType(argv[1], pixelType, componentType, dimensionType);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT(componentType, pixelType, dimensionType, argc, argv);
    }






