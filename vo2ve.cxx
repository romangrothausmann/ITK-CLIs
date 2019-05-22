////program to convert a voxel representation of a simplical 1-complex (e.g. voxel-skeleton consisting only of vertices and lines) into a vector/graph representation
///// ITK only version of vo-skel2poly.cxx (VTK, predecessor), similar but not equivalent 
//01: based on count_neighbours.cxx


/**************************************************************************
NOTE: does not resolve if a branch-point is no longer necessary/apropriate,
e.g. does not resovle this case (neither does prune-ends of skeleton analyzer in fiji):

this:     becomes:   ought to be:
     1
     2
     3         2
 .2233     .2233      .2222
      2         2          2
       2.        2.         2.


this matters for vo2ve but not for the EPC (e.g. calculated with imEuler3D.m)
***************************************************************************/

////ToDo:
// - resovle branchpoints correctly


#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>

///for fill-holes check BEGIN
#include <itkBinaryNotImageFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkShapeOpeningLabelMapFilter.h>//takes a LabelMap as input wheras LabelShapeOpeningImageFilter takes a labeled image as input, LabelImageToLabelMapFilter converts such an image into a LabelMap
#include <itkLabelMapMaskImageFilter.h>
///for fill-holes check END

#include "filter/external/itkCountNeighborsImageFilter/itkCountNeighborsImageFilter.h"
#include <itkChangeLabelImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryImageToStatisticsLabelMapFilter.h>
#include <itkShapeLabelObject.h>
#include <itkLabelMap.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkConstantBoundaryCondition.h>
#include <itkConstNeighborhoodIterator.h>

#include <itkMesh.h>
#include <itkLineCell.h>
#include <itkMeshFileWriter.h>
#include <itkVTKPolyDataMeshIO.h>
#include <itkMetaDataObject.h>



template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef uint8_t  OutputPixelType;
    typedef uint32_t  LabelPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<LabelPixelType, Dimension>  LabelImageType;
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

    const typename InputImageType::Pointer& input= reader->GetOutput();

        {
        ////taken from Modules/Filtering/LabelMap/include/itkBinaryFillholeImageFilter.hxx
        int m_FullyConnected= 0;
        InputPixelType m_ForegroundValue = itk::NumericTraits<InputPixelType>::max();
        InputPixelType backgroundValue = itk::NumericTraits<InputPixelType>::ZeroValue();
        if ( m_ForegroundValue == backgroundValue )
            {
            // current background value is already used for foreground value
            // choose another one
            backgroundValue = itk::NumericTraits<InputPixelType>::max();
            }

        typedef itk::BinaryNotImageFilter< InputImageType > NotType;
        typename NotType::Pointer notInput = NotType::New();
        notInput->SetInput(input);
        notInput->SetForegroundValue( m_ForegroundValue );
        notInput->SetBackgroundValue( backgroundValue );
        notInput->SetReleaseDataFlag( true );
        FilterWatcher watcher0(notInput);

        typedef typename itk::BinaryImageToShapeLabelMapFilter< InputImageType > LabelizerType;
        typename LabelizerType::Pointer labelizer = LabelizerType::New();
        labelizer->SetInput( notInput->GetOutput() );
        labelizer->SetInputForegroundValue( m_ForegroundValue );
        labelizer->SetOutputBackgroundValue( backgroundValue );
        labelizer->SetFullyConnected( m_FullyConnected );
        FilterWatcher watcher1(labelizer);

        typedef typename LabelizerType::OutputImageType                  LabelMapType;
        typedef typename itk::ShapeOpeningLabelMapFilter< LabelMapType > OpeningType;
        typename OpeningType::Pointer opening= OpeningType::New();
        opening->SetInput(labelizer->GetOutput());
        opening->SetAttribute(LabelMapType::LabelObjectType::NUMBER_OF_PIXELS_ON_BORDER);
        opening->SetLambda(1);
        FilterWatcher watcher2(opening);

        // invert the image during the binarization
        typedef typename itk::LabelMapMaskImageFilter< LabelMapType, OutputImageType > BinarizerType;
        typename BinarizerType::Pointer binarizer = BinarizerType::New();
        binarizer->SetInput(opening->GetOutput());
        binarizer->SetLabel(backgroundValue);
        binarizer->SetNegated(true);
        binarizer->SetBackgroundValue(m_ForegroundValue);
        binarizer->SetFeatureImage(input);//will cause another read if reader->ReleaseDataFlagOn();
        FilterWatcher watcher3(binarizer);
        try{
            binarizer->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }

        const typename OutputImageType::Pointer& output= binarizer->GetOutput();

        //// ShapeOpeningLabelMapFilter calls std::map::erase(key) which seems not only to delete the contents but also reduces the total count of elements, so GetNumberOfLabelObjects() reflects the actually contained elements even if not labelled consecutively

        static long int noh= labelizer->GetOutput()->GetNumberOfLabelObjects() - opening->GetOutput()->GetNumberOfLabelObjects();

        if(noh)
            std::cerr << "# of holes removed: " << noh << std::endl << "This generally means that the input does not repesent a simplical 1-complex!!! (Forgot fill_holes-filter before skeletonization?)" << std::endl;
        else
            std::cerr << "No holes found. All fine." << std::endl;
        }

    typename InputImageType::SizeType radius;
    radius.Fill(1);//only regard 26 connectivity

    typedef itk::CountNeighborsImageFilter<InputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);
    filter->SetRadius(radius);
    filter->SetCountNonZero();
    filter->SetValueOfInterest(static_cast<InputPixelType>(atof(argv[3])));
    //filter->ReleaseDataFlagOn();//will be needed multiple times!
    //filter->InPlaceOn();//not available

    FilterWatcher watcher1(filter);
    try{
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    const typename OutputImageType::Pointer& output= filter->GetOutput();


    const OutputPixelType fg= 1;
    const OutputPixelType bg= 0;

    //// get all verts that are neither end-nodes (NN=1) nor primitive connection nodes (NN=2)
    typedef itk::ChangeLabelImageFilter<OutputImageType, OutputImageType> ChangeLabType;
    typename ChangeLabType::Pointer ch= ChangeLabType::New();
    ch->SetInput(filter->GetOutput());
    ch->SetChange(2, 0);//delete connecting nodes of branches

    typedef itk::BinaryThresholdImageFilter<OutputImageType, LabelImageType> ThrFilterType; //thr needs to output LabelImageType such that sd rund with a big cbl is inside value range
    typename ThrFilterType::Pointer thr= ThrFilterType::New();
    thr->SetInput(ch->GetOutput());
    thr->SetLowerThreshold(1);//all non-connecting nodes of any degree
    thr->SetOutsideValue(bg);
    thr->SetInsideValue(fg);
    thr->ReleaseDataFlagOn();
    FilterWatcher watcherThr(thr);

    typedef itk::BinaryImageToStatisticsLabelMapFilter<LabelImageType, OutputImageType> AnaFilterType;
    typename AnaFilterType::Pointer ana= AnaFilterType::New();
    ana->SetInput(thr->GetOutput());
    ana->SetFeatureImage(output);
    ana->SetInputForegroundValue(fg);
    ana->FullyConnectedOn();//expecting skeletons generated by itkBinaryThinningImageFilter3D which are not fully connected so fully connectivity is needed here
    ana->ReleaseDataFlagOn();

    FilterWatcher watcher2(ana);
    try{
        ana->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }
    std::cerr << "# of branching nodes: " << ana->GetOutput()->GetNumberOfLabelObjects() << std::endl;

    //// create mesh to save in a VTK-file
    typedef typename itk::Mesh<float, Dimension>  MeshType;
    typename MeshType::Pointer  mesh = MeshType::New();

    typedef typename AnaFilterType::OutputImageType LabelMapType;
    typedef typename AnaFilterType::OutputImageType::LabelObjectType LabelObjectType;
    typedef typename AnaFilterType::OutputImageType::LabelType LabelType;//default: SizeValueType

    typename LabelMapType::Pointer labelMap = ana->GetOutput();
    const LabelObjectType* labelObject;
    LabelType cbl= labelMap->GetNumberOfLabelObjects() + 1;

    for(LabelType label= 0; label < labelMap->GetNumberOfLabelObjects(); label++){

        labelObject= labelMap->GetNthLabelObject(label);
        typename LabelObjectType::ConstIndexIterator lit(labelObject);

        mesh->SetPoint(label, labelObject->GetCentroid());
        mesh->SetPointData(label, labelObject->GetMean());//make sure filter->ReleaseDataFlagOn() is NOT set!
        }

    typename MeshType::PointIdentifier pointIndex= mesh->GetNumberOfPoints();
    std::cerr << "# of mesh points: " << mesh->GetNumberOfPoints() << std::endl;

    typename LabelImageType::Pointer bpi;
        {//scoped for better consistency
        typedef itk::LabelMapToLabelImageFilter<LabelMapType, LabelImageType> LMtLIType;
        typename LMtLIType::Pointer lmtli = LMtLIType::New();
        lmtli->SetInput(labelMap);
        FilterWatcher watcher3(lmtli);
        try{
            lmtli->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }
        typename LabelImageType::Pointer li= lmtli->GetOutput();
        li->DisconnectPipeline();//essential, otherwise changes to thr will triger a re-update of ana and then lmtli

        //// get connecting nodes
        thr->SetInput(filter->GetOutput());
        thr->SetLowerThreshold(2);//only connecting nodes of branches
        thr->SetUpperThreshold(2);//only connecting nodes of branches
        thr->SetInsideValue(cbl);

        //// add to bpi the connecting nodes with a label of cbl
        typedef itk::AddImageFilter<LabelImageType, LabelImageType, LabelImageType> AddType;
        typename AddType::Pointer adder = AddType::New();
        adder->SetInput1(thr->GetOutput());
        adder->SetInput2(li);
        FilterWatcher watcher4(adder);

        try{
            adder->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }
        bpi= adder->GetOutput();
        bpi->DisconnectPipeline();
        }


    ana->SetInput(bpi);
    ana->SetInputForegroundValue(cbl);

    try{
        ana->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }
    std::cerr << "# of branches: " << ana->GetOutput()->GetNumberOfLabelObjects() << std::endl;

    typedef typename MeshType::CellType CellType;
    typedef typename itk::LineCell<CellType> LineType;
    typedef typename CellType::CellAutoPointer  CellAutoPointer;

    size_t cellId=0;
    typename MeshType::PointType mP;
    for(LabelType label= 0; label < labelMap->GetNumberOfLabelObjects(); label++){

        labelObject= labelMap->GetNthLabelObject(label);
        typename LabelObjectType::ConstIndexIterator lit(labelObject);

        typename LabelImageType::RegionType region= bpi->GetLargestPossibleRegion();
        typedef itk::ConstantBoundaryCondition<LabelImageType> BoundaryConditionType;//zero by default
        typedef itk::ConstNeighborhoodIterator<LabelImageType, BoundaryConditionType> NeighborhoodIteratorType;//ConstantBoundaryCondition essential for not-padded images, default is ZeroFluxNeumannBoundaryCondition
        typename NeighborhoodIteratorType::RadiusType radius;
        radius.Fill(1); //26-connectivity
        NeighborhoodIteratorType tit(radius, bpi, region);

        tit.SetLocation(lit.GetIndex());//lit initial position can be on any pixel of the label

        ////find one branch end
        typename NeighborhoodIteratorType::IndexType lastIndex;
        typename NeighborhoodIteratorType::IndexType lastVIndex;
        lastIndex.Fill(-1);
        lastVIndex.Fill(-1);
        LabelPixelType v= 0;
        bool foundEnd= false;
        while(!foundEnd){
            foundEnd=true;
            typename NeighborhoodIteratorType::OffsetType nextMove;
            nextMove.Fill(0);
            for(unsigned int i = 0; i < tit.Size(); ++i){
                if(i == tit.GetCenterNeighborhoodIndex())
                    continue;//skips center pixel
                if(tit.GetIndex(i) == lastIndex)
                    continue;//skips last visited pixel
                LabelPixelType t= tit.GetPixel(i);
                if(t){
                    if(t == cbl){
                        nextMove= tit.GetOffset(i);
                        foundEnd=false;//this is not an end pixel
                        }
                    else{
                        v= t;
                        lastVIndex= tit.GetIndex(i);
                        }
                    }
                }
            lastIndex= tit.GetIndex();//it.GetIndex() == it.GetIndex(it.GetCenterNeighborhoodIndex())
            tit+= nextMove;
            }


        typename MeshType::PointIdentifier pPointIndex= v-1;//join with bp, bpi value > 0 is index-1 (-1 because labelMap object 0 is mapped to value 1 by LabelMapToLabelImageFilter
        bool foundOtherEnd= false;
        lastIndex= lastVIndex;//for skipping last found bp in first run to find other bp
        while(!foundOtherEnd){

            labelMap->TransformIndexToPhysicalPoint(tit.GetIndex(), mP);
            mesh->SetPoint(pointIndex, mP);
            mesh->SetPointData(pointIndex, labelObject->GetMean());//only connecting nodes of branches

            typename CellType::CellAutoPointer line;
            line.TakeOwnership(new LineType);
            line->SetPointId(0, pPointIndex);
            line->SetPointId(1, pointIndex);
            mesh->SetCell(cellId, line);
            //mesh->SetCellData(cellId, label);//this would be logical but turns out to  be wrong if VTK-file is loaded by e.g. paraview

            pPointIndex= pointIndex;
            pointIndex++;
            cellId++;

            foundOtherEnd=true;
            typename NeighborhoodIteratorType::OffsetType nextMove;
            nextMove.Fill(0);
            for(unsigned int i = 0; i < tit.Size(); ++i){
                if(i == tit.GetCenterNeighborhoodIndex())
                    continue;//skips center pixel
                if(tit.GetIndex(i) == lastIndex)
                    continue;//skips last visited pixel
                LabelPixelType t= tit.GetPixel(i);
                if(t){
                    if(t == cbl){
                        nextMove= tit.GetOffset(i);
                        foundOtherEnd=false;//this is not an end pixel
                        }
                    else
                        v= t;
                    }
                }
            lastIndex= tit.GetIndex();//it.GetIndex() == it.GetIndex(it.GetCenterNeighborhoodIndex())
            tit+= nextMove;
            }

        typename CellType::CellAutoPointer line;
        line.TakeOwnership(new LineType);
        line->SetPointId(0, pPointIndex);
        line->SetPointId(1, v-1);
        mesh->SetCell(cellId, line);
        cellId++;


        mesh->SetCellData(label, label+1);//oddity of ITK-mesh: consecutive line-cells are joined to form a single polyline-cell
        fprintf(stderr, "\r%5.1f%% (%ld)", ((label+1.0) * 100.0 / labelMap->GetNumberOfLabelObjects()), label+1);
        }

    std::cerr << std::endl << "# of mesh points: " << mesh->GetNumberOfPoints() << std::endl;
    std::cerr << "# of mesh cells: " << mesh->GetNumberOfCells() << std::endl;

    int64_t EPC= mesh->GetNumberOfPoints() - mesh->GetNumberOfCells();
    std::cerr << std::endl << "EPC of mesh: " << EPC << std::endl;

    typedef typename itk::MeshFileWriter<MeshType> MeshWriterType;
    typename MeshWriterType::Pointer mwriter = MeshWriterType::New();

    FilterWatcher watcherMO(mwriter);
    mwriter->SetFileName(argv[2]);
    mwriter->SetInput(mesh);

    itk::VTKPolyDataMeshIO::Pointer mio= itk::VTKPolyDataMeshIO::New();
    itk::MetaDataDictionary & metaDic= mio->GetMetaDataDictionary();
    itk::EncapsulateMetaData<std::string>(metaDic, "pointScalarDataName", "Degree");
    itk::EncapsulateMetaData<std::string>(metaDic, "cellScalarDataName", "BranchId");
    mwriter->SetMeshIO(mio);

    try{
        mwriter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
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
    int res= EXIT_FAILURE;
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
    int res= EXIT_FAILURE;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    // case itk::ImageIOBase::CHAR:{         // int8_t
    //     typedef char InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::USHORT:{       // uint16_t
    //     typedef unsigned short InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::SHORT:{        // int16_t
    //     typedef short InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::UINT:{         // uint32_t
    //     typedef unsigned int InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::INT:{          // int32_t
    //     typedef int InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::ULONG:{        // uint64_t
    //     typedef unsigned long InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::LONG:{         // int64_t
    //     typedef long InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::FLOAT:{        // float32
    //     typedef float InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::DOUBLE:{       // float64
    //     typedef double InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
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
                  << " Input_Image"
                  << " Output_VTK-file"
                  << " fg"
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






