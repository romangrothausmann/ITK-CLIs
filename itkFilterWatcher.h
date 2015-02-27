/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFilterWatcher.h,v $
  Language:  C++
  Date:      $Date: 2007-01-29 14:42:11 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkFilterWatcher_h
#define _itkFilterWatcher_h

#include "itkCommand.h"
#include "itkProcessObject.h"
#include <time.h>
#include <stdio.h>

// The following class is a convenience  to watch the progress of a filter

class FilterWatcher
{
 public:
  FilterWatcher(itk::ProcessObject* o, const char *comment="")
    {
      m_Start = 0; m_End = 0; m_Process = o; m_Steps = 0; m_Comment = comment;
      m_TestAbort = false;
      m_Quiet = true;
      m_ReportTime = true;

      itk::SimpleMemberCommand<FilterWatcher>::Pointer startFilterCommand;
      itk::SimpleMemberCommand<FilterWatcher>::Pointer endFilterCommand;
      itk::SimpleMemberCommand<FilterWatcher>::Pointer progressFilterCommand;
      itk::SimpleMemberCommand<FilterWatcher>::Pointer iterationFilterCommand;
      itk::SimpleMemberCommand<FilterWatcher>::Pointer abortFilterCommand;
  
      startFilterCommand =    itk::SimpleMemberCommand<FilterWatcher>::New();
      endFilterCommand =      itk::SimpleMemberCommand<FilterWatcher>::New();
      progressFilterCommand = itk::SimpleMemberCommand<FilterWatcher>::New();
      iterationFilterCommand = itk::SimpleMemberCommand<FilterWatcher>::New();
      abortFilterCommand = itk::SimpleMemberCommand<FilterWatcher>::New();

      startFilterCommand->SetCallbackFunction(this, &FilterWatcher::StartFilter);
      endFilterCommand->SetCallbackFunction(this, &FilterWatcher::EndFilter);
      progressFilterCommand->SetCallbackFunction(this, &FilterWatcher::ShowProgress);
      iterationFilterCommand->SetCallbackFunction(this, &FilterWatcher::ShowIteration);
      abortFilterCommand->SetCallbackFunction(this, &FilterWatcher::ShowAbort);
      m_Process->AddObserver(itk::StartEvent(), startFilterCommand);
      m_Process->AddObserver(itk::EndEvent(), endFilterCommand);
      m_Process->AddObserver(itk::ProgressEvent(), progressFilterCommand);
      m_Process->AddObserver(itk::IterationEvent(), iterationFilterCommand);
      m_Process->AddObserver(itk::AbortEvent(), abortFilterCommand);
    }

  virtual ~FilterWatcher() {}

  virtual void ShowProgress(){
    m_Steps++;
    //if (!m_Quiet)
    {
      fprintf(stderr, "\r%s progress: %5.1f%%", m_Process->GetNameOfClass(), 100.0 * m_Process->GetProgress());
      std::cerr.flush();
    }
    if (m_TestAbort)
      {
	if (m_Process->GetProgress() > .03)
	  {
	    m_Process->AbortGenerateDataOn();
	  }
      }
  }

  virtual void ShowAbort(){
    std::cerr << std::endl << "      ABORT" << std::endl << std::flush;
  }

  virtual void ShowIteration(){
    std::cerr << " # " << std::flush;
    m_Iterations++;
  }

  virtual void StartFilter(){
    m_Steps = 0;
    m_Iterations = 0;
    m_Start = ::clock();
    m_tStart = ::time(NULL);
    if (!m_Quiet){
      std::cerr << "-------- Start " << m_Process->GetNameOfClass()
		<< " \"" << m_Comment << "\" "
		<< m_Process
		<< (m_Quiet ? "Progress Quiet " : "Progress ")
		<< std::flush;
    }
    else{
      //std::cerr << "Executing " << m_Process->GetNameOfClass() << " " << std::flush;
      std::cerr << m_Process->GetNameOfClass() << " " << std::flush;
    }
  }


  virtual void EndFilter(){
    std::cerr << " done.";
    if (m_ReportTime){
        m_End = ::clock();
        m_tEnd = ::time(NULL);
 	//std::cerr << std::endl << std::flush;
	fprintf(stderr, " (real= %.1fs, CPUs= %.1fs)", ::difftime(m_tEnd, m_tStart), static_cast<double>(m_End - m_Start) / CLOCKS_PER_SEC);
	/* std::cerr << " (real= " */
        /*           << ::difftime(m_tEnd, m_tStart) << " s, sys= " */
        /*           << static_cast<double>(m_End - m_Start) / CLOCKS_PER_SEC */
        /*           << " s)."; */
        /* std::cerr << std::endl << std::endl */
        /*           << "-------- End " << m_Process->GetNameOfClass() */
        /*           << " \"" << m_Comment << "\" " */
        /*           << m_Process << std::flush; */
      }
    ////if there is no progess info, just do nothing!
    /* if (m_Steps < 1){ */
    /*   /\* if (!m_Quiet)  *\/ */
    /*   /\* 	std::cerr << m_Process->GetNameOfClass() << " has no progress info." << std::endl << std::flush; *\/ */
    /*   /\* itkExceptionMacro ("Filter does not have progress."); //returns immediately! *\/ */
    /*   //do not throw an execption if there is no progress! */
    /*   std::cerr << " has no progress info." << std::endl << std::flush; */
    /* } */
    std::cerr << std::endl << std::flush;   
  }
  

  const char *GetNameOfClass () {return m_Process->GetNameOfClass();}

  void QuietOn() {m_Quiet = true;};
  void QuietOff() {m_Quiet = false;};
  void TestAbortOn() {m_TestAbort = true;};
  void TestAbortOff() {m_TestAbort = false;};
  void ReportTimeOn() {m_ReportTime= true;};
  void ReportTimeOff() {m_ReportTime= false;};

 protected:
  clock_t m_Start;
  clock_t m_End;
  time_t m_tStart;
  time_t m_tEnd;
  int m_Steps;
  int m_Iterations;
  bool m_Quiet;
  bool m_TestAbort;
  bool m_ReportTime;
  std::string m_Comment;
  itk::ProcessObject::Pointer m_Process;
 private:
  FilterWatcher(); // Purposely not implemented
};

#endif
