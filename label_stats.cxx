////program for itkLabelStatisticsImageFilter
//01: based on template_2inputs.cxx and analyse_labels.cxx


#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkLabelImageToStatisticsLabelMapFilter.h>
#include <itkStatisticsLabelObject.h>
#include <itkLabelMap.h>



template<typename InputComponent1, typename TypeInputComponentType2, typename InputPixelType1, typename InputPixelType2, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef itk::Image<InputPixelType1, Dimension>  InputImageType1;
    typedef itk::Image<InputPixelType2, Dimension>  InputImageType2;


    typedef itk::ImageFileReader<InputImageType1> ReaderType1;
    typename ReaderType1::Pointer reader1 = ReaderType1::New();

    reader1->SetFileName(argv[1]);
    reader1->ReleaseDataFlagOn();
    FilterWatcher watcherI1(reader1);
    watcherI1.QuietOn();
    watcherI1.ReportTimeOn();
    try{
        reader1->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    const typename InputImageType1::Pointer& input1= reader1->GetOutput();


    typedef itk::ImageFileReader<InputImageType2> ReaderType2;
    typename ReaderType2::Pointer reader2 = ReaderType2::New();

    reader2->SetFileName(argv[2]);
    reader2->ReleaseDataFlagOn();
    FilterWatcher watcherI2(reader2);
    watcherI2.QuietOn();
    watcherI2.ReportTimeOn();
    try{
        reader2->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    const typename InputImageType2::Pointer& input2= reader2->GetOutput();


    typedef itk::LabelImageToStatisticsLabelMapFilter<InputImageType1, InputImageType2> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput1(reader1->GetOutput());
    filter->SetInput2(reader2->GetOutput()); //SetFeatureImage
    filter->ReleaseDataFlagOn();
    bool cp= atoi(argv[3]);
    filter->SetComputePerimeter(cp);
    filter->SetComputeFeretDiameter(0);
    filter->SetComputeHistogram(0);
    // filter->SetNumberOfBins();

    FilterWatcher watcher1(filter);
    try{
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    typedef typename FilterType::OutputImageType LabelMapType;
    typedef typename FilterType::OutputImageType::LabelObjectType LabelObjectType;
    typedef typename FilterType::OutputImageType::LabelType LabelType;//default: InputImageType::PixelType

    // then we can read the attribute values we're interested in. The BinaryImageToShapeLabelMapFilter
    // produce consecutive labels, so we can use a for loop and GetLabelObject() method to retrieve
    // the label objects. If the labels are not consecutive, the GetNthLabelObject() method must be
    // use instead of GetLabelObject(), or an iterator on the label object container of the label map.
    //// list of quantities see: http://www.itk.org/Doxygen/html/classitk_1_1ShapeLabelObject.html
    typename LabelMapType::Pointer labelMap = filter->GetOutput();
    std::cout << "#index";
    for (unsigned int i= 0; i < Dimension; i++)
        std::cout << "\tbbIndex_" << i+1;
    for (unsigned int i= 0; i < Dimension; i++)
        std::cout << "\tbbSize_" << i+1;
    for (unsigned int i= 0; i < Dimension; i++)
        std::cout << "\tbary_" << i+1;
    for (unsigned int i= 0; i < Dimension; i++)
        std::cout << "\tellALength_" << i+1;
    for (unsigned int i= 0; i < Dimension; i++)
        for (unsigned int j= 0; j < Dimension; j++)
            std::cout << "\tellAOri_" << i+1 << j+1;
    std::cout << "\tvoxel\tVphy";
    if (cp)
        std::cout << "\tAphy";

    std::cout << "\tmin\tmax\tmean\tstd\tvar\tsum";
    
    for (unsigned int i= 0; i < Dimension; i++)
        std::cout << "\tgrav_" << i+1;
    for (unsigned int i= 0; i < Dimension; i++)
        std::cout << "\tmin_" << i+1;
    for (unsigned int i= 0; i < Dimension; i++)
        std::cout << "\tmax_" << i+1;

    std::cout << std::endl;

    const LabelObjectType* labelObject;
    for(LabelType label= 0; label < labelMap->GetNumberOfLabelObjects(); label++){//SizeValueType == LabelType //GetNthLabelObject starts with 0 and ends at GetNumberOfLabelObjects()-1!!!

        labelObject= labelMap->GetNthLabelObject(label);//GetNthLabelObject essential in case of not consecutive labels! (compare: http://www.itk.org/Doxygen47/html/classitk_1_1BinaryImageToShapeLabelMapFilter.html and http://www.itk.org/Doxygen47/html/classitk_1_1LabelMap.html)
    	// check for new measures (e.g. OBB): https://itk.org/Doxygen/html/classitk_1_1StatisticsLabelObject.html
        std::cout
            << +labelObject->GetLabel() << "\t";//essential in case of not consecutive labels!
        for (unsigned int i= 0; i < Dimension; i++)
            std::cout << +labelObject->GetBoundingBox().GetIndex()[i] << "\t";
        for (unsigned int i= 0; i < Dimension; i++)
            std::cout << +labelObject->GetBoundingBox().GetSize()[i] << "\t";
        for (unsigned int i= 0; i < Dimension; i++)
            std::cout << +labelObject->GetCentroid()[i] << "\t";
        for (unsigned int i= 0; i < Dimension; i++)
            std::cout << +labelObject->GetEquivalentEllipsoidDiameter()[i] << "\t";
        for (unsigned int i= 0; i < labelObject->GetPrincipalAxes().GetVnlMatrix().columns(); i++)
            for (unsigned int j= 0; j < labelObject->GetPrincipalAxes().GetVnlMatrix().rows(); j++)
                std::cout << +labelObject->GetPrincipalAxes()[i][j] << "\t";
        std::cout
            //<< labelObject->GetFeretDiameter() << "\t"
            //<< labelObject->GetElongation() << "\t"
            //<< labelObject->GetRoundness() << "\t" // Lehmann2008 Sec. 8.2: Roundness here is calc. like Sphericity: Psi= Ss(V)/S
            //<< labelObject->GetFlatness() << "\t"
            //<< labelObject->GetEquivalentSphericalRadius() << "\t"
            //<< labelObject->GetEquivalentSphericalPerimeter() << "\t"
            << +labelObject->GetNumberOfPixels() << "\t"
            << +labelObject->GetPhysicalSize();
        if (cp)
            std::cout << "\t" << +labelObject->GetPerimeter();//takes physical voxel size into account!!!
	
	std::cout
	    << "\t" << +labelObject->GetMinimum()
	    << "\t" << +labelObject->GetMaximum()
	    << "\t" << +labelObject->GetMean()
	    << "\t" << +labelObject->GetStandardDeviation()
	    << "\t" << +labelObject->GetVariance()
	    // << "\t" << +labelObject->GetKurtosis()
	    // << "\t" << +labelObject->GetSkewness()
	    // << "\t" << +labelObject->GetWeightedElongation()
	    // << "\t" << +labelObject->GetWeightedFlatness()
	    // << "\t" << +labelObject->GetHistogram()
	    // << "\t" << +labelObject->GetMedian() // needs histogram
	    << "\t" << +labelObject->GetSum();

        for (unsigned int i= 0; i < Dimension; i++)
            std::cout << "\t" << +labelObject->GetCenterOfGravity()[i];
        for (unsigned int i= 0; i < Dimension; i++)
            std::cout << "\t" << +labelObject->GetMinimumIndex()[i];
        for (unsigned int i= 0; i < Dimension; i++)
            std::cout << "\t" << +labelObject->GetMaximumIndex()[i];
	
        std::cout << std::endl;
        }

    return EXIT_SUCCESS;

    }


template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1, typename InputPixelType2>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    // case 1:
    //     res= DoIt<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2, 1>(argc, argv);
    //     break;
    case 2:
        res= DoIt<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2, 2>(argc, argv);
        break;
    case 3:
        res= DoIt<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2, 3>(argc, argv);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1>
int dispatch_pT2(itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType2){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType2 InputPixelType2;
        res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
        } break;
    // case itk::ImageIOBase::RGB:{
    //     typedef itk::RGBPixel<InputComponentType2> InputPixelType2;
    //     res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::RGBA:{
    //     typedef itk::RGBAPixel<InputComponentType2> InputPixelType2;
    //     res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::COMPLEX:{
    //     typedef std::complex<InputComponentType2> InputPixelType2;
    //     res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::VECTOR:{
    //     typedef itk::VariableLengthVector<InputComponentType2> InputPixelType2;
    //     res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
    //     } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType1, typename InputComponentType2>
int dispatch_pT1(itk::ImageIOBase::IOPixelType pixelType1, itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType1){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType1 InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType1>
int dispatch_cT2(itk::ImageIOBase::IOComponentType componentType2, itk::ImageIOBase::IOPixelType pixelType1, itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType2){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::FLOAT:{        // float32
        typedef float InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::DOUBLE:{       // float64
        typedef double InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    default:
        std::cerr << "unknown component type" << std::endl;
        break;
        }//switch
    return res;
    }

int dispatch_cT1(itk::ImageIOBase::IOComponentType componentType1, itk::ImageIOBase::IOComponentType componentType2, itk::ImageIOBase::IOPixelType pixelType1, itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType1){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
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

    if(!imageReader->GetImageIO()->CanStreamRead())
        std::cerr << "Cannot stream the reading of the input. Streaming will be inefficient!" << std::endl;

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
    if ( argc != 4 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Label_Image"
                  << " Grey_Image"
                  << " <bool: calculate perimeter>"
                  << std::endl;

        return EXIT_FAILURE;
        }

    itk::ImageIOBase::IOPixelType pixelType1;
    itk::ImageIOBase::IOPixelType pixelType2;
    typename itk::ImageIOBase::IOComponentType componentType1;
    typename itk::ImageIOBase::IOComponentType componentType2;
    size_t dimensionType1;
    size_t dimensionType2;


    try {
        GetImageType(argv[1], pixelType1, componentType1, dimensionType1);
        GetImageType(argv[2], pixelType2, componentType2, dimensionType2);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    if (dimensionType1 != dimensionType2){
        std::cout << "Input images need to be of the same dimension." << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT1(componentType1, componentType2, pixelType1, pixelType2, dimensionType1, argc, argv);
    }






