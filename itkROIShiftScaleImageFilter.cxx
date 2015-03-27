#ifndef __itkROIShiftScaleImageFilter_hxx
#define __itkROIShiftScaleImageFilter_hxx
 
#include "itkROIShiftScaleImageFilter.h"

namespace itk{

    template<typename TInputImage, typename TOutputImage>
    ROIShiftScaleImageFilter<TInputImage, TOutputImage>
    ::ROIShiftScaleImageFilter(){
	}

    // template<typename TInputImage, typename TOutputImage>
    // void
    // ROIShiftScaleImageFilter<TInputImage, TOutputImage>
    // ::SetROI(InputImageRegionType extractRegion)
    // 	{
    // 	m_ROI = extractRegion;
    // 	this->Modified();
    // 	}

    template<typename TInputImage, typename TOutputImage>
    void ROIShiftScaleImageFilter<TInputImage, TOutputImage>::GenerateData(){

	const typename TInputImage::ConstPointer input = this->GetInput();
	const typename TInputImage::RegionType  region = input->GetRequestedRegion();
	const typename TInputImage::SizeType      size = region.GetSize();
 
	typename TOutputImage::Pointer output = this->GetOutput();
 

	typename ExtractFilterType::Pointer extract= ExtractFilterType::New();
	extract->SetInput(input);
	extract->SetExtractionRegion(m_ROI);
	extract->SetDirectionCollapseToIdentity(); // This is required.

	try{ 
	    extract->Update();
	    }
	catch(itk::ExceptionObject &ex){ 
	    std::cerr << ex << std::endl;
	    }

	typename StatImageFilterType::Pointer stat= StatImageFilterType::New();
	stat->SetInput(extract->GetOutput());
	//stat->GetOutput()->SetRequestedRegion(ROI); //does not work for StatisticsImageFilter
	//stat->GenerateInputRequestedRegion(); //does not help either
	try { 
	    stat->Update();
	    }
	catch(itk::ExceptionObject &ex){ 
	    std::cerr << ex << std::endl;
	    }

	//std::cerr << "Min: " << itk::ImageIOBase::IOPixelType(stat->GetMinimum()) << " Max: " << itk::ImageIOBase::IOPixelType(stat->GetMaximum()) << " Mean: " << stat->GetMean() << " Std: " << stat->GetSigma() << " Variance: " << stat->GetVariance() << " Sum: " << stat->GetSum() << std::endl; //stat->GetMaximum() returns PixelType which needs to be casted for proper output, using itk::ImageIOBase::IOPixelType here somehow seems to work right...

	typename SSImageFilterType::Pointer ss= SSImageFilterType::New();
	ss->SetShift(static_cast<typename SSImageFilterType::RealType>(m_DesiredMean - stat->GetMean() * m_DesiredStd/stat->GetSigma()));
	ss->SetScale(static_cast<typename SSImageFilterType::RealType>(m_DesiredStd/stat->GetSigma()));
	//std::cerr << " Shift: " << ss->GetShift() << " Scale: " << ss->GetScale() << std::endl;
	ss->SetInput(input);

	try { 
	    ss->Update();
	    }
	catch(itk::ExceptionObject &ex){ 
	    std::cerr << ex << std::endl;
	    }

	output->Graft(ss->GetOutput());

	}
 
    template<typename TInputImage, typename TOutputImage>
    void
    ROIShiftScaleImageFilter<TInputImage, TOutputImage>
    ::PrintSelf(std::ostream & os, Indent indent) const
	{
	Superclass::PrintSelf(os, indent);

	os << indent << "ROI: " << m_ROI << std::endl;
	os << indent << "DesiredMean: " << m_DesiredMean << std::endl;
	os << indent << "DesiredStd: " << m_DesiredStd << std::endl;
	}


    }// end namespace
 
#endif //__itkROIShiftScaleImageFilter_hxx
