////program for itkSpeedFunctionToPathFilter (Fast Marching Minimal Path Extraction)
//01: based on template.cxx and vo-img2v-skel_01.cxx


#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkSpeedFunctionToPathFilter.h>
#include <itkGradientDescentOptimizer.h>
#include <itkPathIterator.h>
#include <itkPolyLineParametricPath.h>
#include <itkImageFileWriter.h>



template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    const char offset= 4;
    if( argc != offset + 3*Dimension){
        fprintf(stderr, "%d + 3*Dimension = %d parameters are needed!\n", offset-1, offset-1 + 3*Dimension);
        return EXIT_FAILURE;
        }

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
    optimizer->SetNumberOfIterations(1000);

    // Create path filter
    typename PathFilterType::Pointer pathFilter = PathFilterType::New();
    pathFilter->SetInput(input); //needs the image values to be scaled to [0; 1]
    pathFilter->SetCostFunction(cost);
    pathFilter->SetOptimizer(optimizer);
    pathFilter->SetTerminationValue(2.0);

    // Setup path points
    typename PathFilterType::PointType start, end, way;

    for(int j= 0; j < Dimension; j++){
        start[j]= atof(argv[j+offset]);
        }

    for(int j= 0; j < Dimension; j++){
        end[j]= atof(argv[j+Dimension+offset]);
        }

    for(int j= 0; j < Dimension; j++){
        way[j]= atof(argv[j+2*Dimension+offset]);
        }

    // Add path information
    typename PathFilterType::PathInfo info;
    info.SetStartPoint(start);
    info.SetEndPoint(end);
    info.AddWayPoint(way);
    pathFilter->AddPathInfo(info);


    // Compute the path
    FilterWatcher watcher(pathFilter); //filter reports no progress so far
    try {
        pathFilter->Update();
        }
    catch (itk::ExceptionObject &ex){
        std::cout << ex << std::endl;
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
    for (unsigned int i=0; i<pathFilter->GetNumberOfOutputs(); i++){

        // Get the path
        typename PathType::Pointer path = pathFilter->GetOutput(i);

        // Check path is valid
        if (path->GetVertexList()->Size() == 0){
            std::cout << "WARNING: Path " << (i+1) << " contains no points!" << std::endl;
            continue;
            }

        // Iterate path and convert to image
        PathIteratorType it(output, path);
        for (it.GoToBegin(); !it.IsAtEnd(); ++it){
            it.Set(itk::NumericTraits<OutputPixelType>::max());
            }
        }


    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
    writer->SetInput(output);
    writer->SetUseCompression(atoi(argv[3]));
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

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
    if ( argc < 4 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress"
                  << " start-point..."
                  << " end-point..."
                  << " way-point..."
                  << std::endl;
        std::cerr << " Point coordinates are expected in physical space!" << std::endl;

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






