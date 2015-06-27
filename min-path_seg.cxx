////program to segment a binray image along a minimal path
//01: based on min-path_fm.cxx and skeletonize_3D.cxx

////ToDo:


#include <string>
#include <sstream>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkParabolicOpenImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkSpeedFunctionToPathFilter.h>
#include <itkMorphologicalDistanceTransformImageFilter.h>
#include <itkGradientDescentOptimizer.h>
#include <itkPathIterator.h>
#include <itkPolyLineParametricPath.h>
#include <itkParabolicDilateImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageFileWriter.h>

#include <itkMesh.h>
#include <itkLineCell.h>
#include <itkMeshFileWriter.h>


template<typename OptimizerType>
void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

    const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

    if(itk::ProgressEvent().CheckEvent(&event))
        fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
    else if(itk::IterationEvent().CheckEvent(&event)){
        OptimizerType* optimizer= dynamic_cast<OptimizerType *>(caller);
        fprintf(stderr, "\r%5d: %7.3f ", optimizer->GetCurrentIteration(), optimizer->GetValue());
        std::cerr << optimizer->GetCurrentPosition();
        }
    else if(itk::EndEvent().CheckEvent(&event))
        std::cerr << std::endl;
    }


template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    const char offset= 7;
    if((argc - offset) % Dimension){
        fprintf(stderr, "%d + n*Dimension  parameters are needed!\n", offset-1);
        return EXIT_FAILURE;
        }

    char* outPrefix= argv[2];
    std::stringstream sss;

    typedef double   SpeedPixelType;
    typedef double   DMPixelType;
    typedef DMPixelType  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<SpeedPixelType, Dimension>  SpeedImageType;
    typedef itk::Image<DMPixelType, Dimension>  DMImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    //reader->ReleaseDataFlagOn();
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
    typedef typename itk::RescaleIntensityImageFilter<InputImageType, InputImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(reader->GetOutput());
    rescaleFilter->SetOutputMinimum(0.0);
    rescaleFilter->SetOutputMaximum(1.0);
    FilterWatcher watcherR(rescaleFilter);

    typedef itk::ParabolicOpenImageFilter<InputImageType, SpeedImageType> SmoothFilterType;
    typename SmoothFilterType::Pointer smoother= SmoothFilterType::New();
    smoother->SetInput(rescaleFilter->GetOutput());
    smoother->SetScale(atof(argv[4]));
    FilterWatcher watcherSM(smoother);


    try {
        smoother->Update();
        }
    catch (itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    const typename SpeedImageType::Pointer& speed= smoother->GetOutput();

        {//// for debugging
        typedef itk::ImageFileWriter<SpeedImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();

        FilterWatcher watcherO(writer);
        sss.str(""); sss << outPrefix << "_speed.mha";
        writer->SetFileName(sss.str().c_str());
        writer->SetInput(speed);
        writer->SetUseCompression(atoi(argv[3]));
        try{
            writer->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }
        }

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
    optimizer->SetNumberOfIterations(atoi(argv[5]));
    optimizer->SetLearningRate(atof(argv[6]));

    itk::CStyleCommand::Pointer eventCallbackITK;
    eventCallbackITK = itk::CStyleCommand::New();
    eventCallbackITK->SetCallback(FilterEventHandlerITK<OptimizerType>);
    optimizer->AddObserver(itk::AnyEvent(), eventCallbackITK);

    // Create path filter
    typename PathFilterType::Pointer pathFilter = PathFilterType::New();
    pathFilter->SetInput(speed); //needs the image values to be scaled to [0; 1]
    pathFilter->SetCostFunction(cost);
    pathFilter->SetOptimizer(optimizer);
    pathFilter->SetTerminationValue(2.0);

    // Setup path points
    typename SpeedImageType::IndexType start, end;
    typename PathFilterType::PointType startP, endP;
    typename PathFilterType::PathInfo info;

    for(int j= 0; j < Dimension; j++){
        start[j]= atoi(argv[j+offset]) - 1;
        }

    for(int j= 0; j < Dimension; j++){
        end[j]= atoi(argv[j+Dimension+offset]) - 1;
        }

    speed->TransformIndexToPhysicalPoint(start, startP);
    speed->TransformIndexToPhysicalPoint(end, endP);
    info.SetStartPoint(startP);
    info.SetEndPoint(endP);
    std::cerr << "S: " << start << " physical coords: "<< startP << std::endl;	
    std::cerr << "E: " << end << " physical coords: "<< endP << std::endl;	

    for(int i= offset + 2*Dimension; i < argc; i+= Dimension){
	typename SpeedImageType::IndexType way;
	typename PathFilterType::PointType wayP;

	for(int j= 0; j < Dimension; j++){
	    way[j]= atoi(argv[i+j]) - 1;
	    }

	speed->TransformIndexToPhysicalPoint(way, wayP);
	info.AddWayPoint(wayP);
	std::cerr << "W: " << way << " physical coords: "<< wayP << std::endl;	
	}

    pathFilter->AddPathInfo(info);
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
    output->SetRegions(speed->GetLargestPossibleRegion());
    output->SetSpacing(speed->GetSpacing());
    output->SetOrigin(speed->GetOrigin());
    output->Allocate();
    output->FillBuffer(itk::NumericTraits<OutputPixelType>::Zero);


    typedef itk::MorphologicalDistanceTransformImageFilter<InputImageType, DMImageType> DMFilterType;
    typename DMFilterType::Pointer dm= DMFilterType::New();
    dm->SetInput(reader->GetOutput());
    //dm->ReleaseDataFlagOn();
    dm->SqrDistOn();
    dm->SetUseImageSpacing(false);

    FilterWatcher watcherDM(dm);
    try{ 
        dm->Update();
        }
    catch(itk::ExceptionObject &ex){ 
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

        {//// for debugging
        typedef itk::ImageFileWriter<DMImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();

        FilterWatcher watcherO(writer);
        sss.str(""); sss << outPrefix << "_dm.mha";
        writer->SetFileName(sss.str().c_str());
        writer->SetInput(dm->GetOutput());
        writer->SetUseCompression(atoi(argv[3]));
        try{
            writer->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }
        }

    // Rasterize path
    typedef itk::PathIterator<OutputImageType, PathType> PathIteratorType;
    typedef itk::PathConstIterator<DMImageType, PathType> PathConstIteratorType;
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
        std::cerr << path->EndOfInput() << std::endl;

        // Iterate path and convert to image
        PathIteratorType it(output, path);
        PathConstIteratorType cit(dm->GetOutput(), path);
        int count= 0;
        for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++cit){//possibly more iterations than path->GetVertexList()->Size() !! can take long with large radii
            it.Set(cit.Get());//mark center with squared hight, later used by ParabolicDilate
            }
        }

    typedef itk::ParabolicDilateImageFilter<DMImageType, DMImageType> PDFilterType;
    typename PDFilterType::Pointer pd= PDFilterType::New();
    pd->SetInput(output);
    pd->SetScale(.5);//0.5: unscaled
    FilterWatcher watcherPD(pd);

    typedef itk::BinaryThresholdImageFilter<DMImageType, OutputImageType> ThrFilterType;
    typename ThrFilterType::Pointer thr= ThrFilterType::New();
    thr->SetInput(pd->GetOutput());
    thr->SetUpperThreshold(0);
    thr->SetOutsideValue(1);
    thr->SetInsideValue(0);
    FilterWatcher watcherThr(thr);

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    sss.str(""); sss << outPrefix << ".mha";
    writer->SetFileName(sss.str().c_str());
    writer->SetInput(thr->GetOutput());
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

    typename MeshType::PointType mP;
    for (unsigned int i=0; i < pathFilter->GetNumberOfOutputs(); i++){

        // Get the path, coords are stored as continous index
        typename PathType::Pointer path = pathFilter->GetOutput(i);
        const typename PathType::VertexListType *vertexList = path->GetVertexList();

        for(unsigned int k = 0; k < vertexList->Size(); k++){
            speed->TransformContinuousIndexToPhysicalPoint(vertexList->GetElement(k), mP);
            mesh->SetPoint(k, mP);
            }
        }

    std::cout << "# of mesh points: " << mesh->GetNumberOfPoints() << std::endl;

    //// create connecting lines
    typedef typename MeshType::CellType CellType;
    typedef typename itk::LineCell<CellType> LineType;
    typedef typename CellType::CellAutoPointer  CellAutoPointer;

    //// from: http://www.itk.org/Doxygen/html/Examples_2DataRepresentation_2Mesh_2Mesh3_8cxx-example.html
    if(mesh->GetNumberOfPoints() > 1){
        const unsigned int numberOfCells = mesh->GetNumberOfPoints() - 1;
        typename CellType::CellAutoPointer line;
        for(size_t cellId=0; cellId < numberOfCells; cellId++){
            line.TakeOwnership(new LineType);
            line->SetPointId(0, cellId);
            line->SetPointId(1, cellId+1);
            mesh->SetCell(cellId, line);
            }
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
    if ( argc < 7 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image_Base"
                  << " compress"
                  << " sigma"
                  << " iterations"
                  << " step-scale"
                  << " start-point..."
                  << " end-point..."
                  << " way-point..."
                  << std::endl;
        std::cerr << " Point coordinates are expected in voxel units (starting with 1)!" << std::endl;
        std::cerr << " step-scale will correspond to distance between points for a speed function ~1 along path" << std::endl;

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






