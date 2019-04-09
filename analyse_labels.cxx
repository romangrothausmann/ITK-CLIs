////program to analyse the labels of a label image
//01: based on analyse_binary.cxx


#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkShapeLabelObject.h>
#include <itkLabelMap.h>




template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;

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

    const typename InputImageType::Pointer& input= reader->GetOutput();


    typedef itk::LabelImageToShapeLabelMapFilter<InputImageType> FilterType;//default output is a LabelMap instanciated with InputImageType::PixelType
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);
    filter->ReleaseDataFlagOn();
    bool cp= atoi(argv[2]);
    filter->SetComputePerimeter(cp);

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
    std::cout << std::endl;

    const LabelObjectType* labelObject;
    for(LabelType label= 0; label < labelMap->GetNumberOfLabelObjects(); label++){//SizeValueType == LabelType //GetNthLabelObject starts with 0 and ends at GetNumberOfLabelObjects()-1!!!

        labelObject= labelMap->GetNthLabelObject(label);//GetNthLabelObject essential in case of not consecutive labels! (compare: http://www.itk.org/Doxygen47/html/classitk_1_1BinaryImageToShapeLabelMapFilter.html and http://www.itk.org/Doxygen47/html/classitk_1_1LabelMap.html)
	// check for new measures (e.g. OBB): https://itk.org/Doxygen/html/classitk_1_1ShapeLabelObject.html
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
        std::cout << std::endl;
        }


    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= 0;
    switch (dimensionType){
    // case 1:
    //     res= DoIt<InputComponentType, InputPixelType, 1>(argc, argv);
    //     break;
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
    if ( argc != 3 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Label_Image"
                  << " <bool: calculate perimeter>"
                  << std::endl;

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






