#ifndef __itkLabelShiftImageFilter_hxx
#define __itkLabelShiftImageFilter_hxx
 
#include "itkLabelShiftImageFilter.h"

namespace itk{

    template<typename TInputImage, typename TOutputImage>
    LabelShiftImageFilter<TInputImage, TOutputImage>
    ::LabelShiftImageFilter(){

	m_AddImageFilter= AddFilterType::New();

	}

    template<typename TInputImage, typename TOutputImage>
    void LabelShiftImageFilter<TInputImage, TOutputImage>
    //::GenerateData(){
    ::ThreadedGenerateData(const OutputImageRegionType & region, ThreadIdType threadId){
	//std::cerr << "Thread " << threadId << " given region: " << region << std::endl;


	const typename TInputImage::ConstPointer input = this->GetInput();
	//const typename TInputImage::RegionType  region = input->GetRequestedRegion();
	//const typename TInputImage::SizeType      size = region.GetSize();
 
	typename TOutputImage::Pointer output = this->GetOutput();
 

	m_AddImageFilter->SetInput1(input);
	m_AddImageFilter->SetConstant2(m_LastMax);
	m_AddImageFilter->GetOutput()->SetRequestedRegion(region);

	m_AddImageFilter->Update();
	output->Graft(m_AddImageFilter->GetOutput()); 

	m_LastMax++;
	}
 
    template<typename TInputImage, typename TOutputImage>
    void
    LabelShiftImageFilter<TInputImage, TOutputImage>
    ::PrintSelf(std::ostream & os, Indent indent) const
	{
	Superclass::PrintSelf(os, indent);

	os << indent << "LastMax: " << m_LastMax << std::endl;
	}


    }// end namespace
 
#endif //__itkLabelShiftImageFilter_hxx
