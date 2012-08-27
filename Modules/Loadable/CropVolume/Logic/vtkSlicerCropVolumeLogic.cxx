/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkSlicerCropVolumeLogic.cxx,v $
  Date:      $Date: 2006/01/06 17:56:48 $
  Version:   $Revision: 1.58 $

=========================================================================auto=*/

// CLI invocation
#include <qSlicerCLIModule.h>
#include <vtkSlicerCLIModuleLogic.h>

// CropLogic includes
#include "vtkSlicerCLIModuleLogic.h"
#include "vtkSlicerCropVolumeLogic.h"
#include "vtkSlicerVolumesLogic.h"

// CropMRML includes

// MRML includes
#include <vtkMRMLAnnotationROINode.h>
#include <vtkMRMLCropVolumeParametersNode.h>
#include <vtkMRMLDiffusionTensorVolumeNode.h>
#include <vtkMRMLDiffusionWeightedVolumeNode.h>
#include <vtkMRMLVectorVolumeNode.h>
#include <vtkMRMLLinearTransformNode.h>
#include <vtkMRMLVolumeNode.h>
#include <vtkMRMLAnnotationROINode.h>

// VTK includes
#include <vtkImageData.h>
#include <vtkImageClip.h>
#include <vtkNew.h>
#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>

// STD includes
#include <cassert>
#include <iostream>

//----------------------------------------------------------------------------
class vtkSlicerCropVolumeLogic::vtkInternal
{
public:
  vtkInternal();

  vtkSlicerVolumesLogic* VolumesLogic;
  vtkSlicerCLIModuleLogic* ResampleLogic;
};

//----------------------------------------------------------------------------
vtkSlicerCropVolumeLogic::vtkInternal::vtkInternal()
{
  this->VolumesLogic = 0;
  this->ResampleLogic = 0;
}

//----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkSlicerCropVolumeLogic, "$Revision: 1.9.12.1 $");
vtkStandardNewMacro(vtkSlicerCropVolumeLogic);

//----------------------------------------------------------------------------
vtkSlicerCropVolumeLogic::vtkSlicerCropVolumeLogic()
{
  this->Internal = new vtkInternal;
}

//----------------------------------------------------------------------------
vtkSlicerCropVolumeLogic::~vtkSlicerCropVolumeLogic()
{
  delete this->Internal;
}

//----------------------------------------------------------------------------
void vtkSlicerCropVolumeLogic::SetVolumesLogic(vtkSlicerVolumesLogic* logic)
{
  this->Internal->VolumesLogic = logic;
}

//----------------------------------------------------------------------------
vtkSlicerVolumesLogic* vtkSlicerCropVolumeLogic::GetVolumesLogic()
{
  return this->Internal->VolumesLogic;
}

//----------------------------------------------------------------------------
void vtkSlicerCropVolumeLogic::SetResampleLogic(vtkSlicerCLIModuleLogic* logic)
{
  this->Internal->ResampleLogic = logic;
}

//----------------------------------------------------------------------------
vtkSlicerCLIModuleLogic* vtkSlicerCropVolumeLogic::GetResampleLogic()
{
  return this->Internal->ResampleLogic;
}

//----------------------------------------------------------------------------
void vtkSlicerCropVolumeLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->vtkObject::PrintSelf(os, indent);
  os << indent << "vtkSlicerCropVolumeLogic:             " << this->GetClassName() << "\n";
}

//----------------------------------------------------------------------------
int vtkSlicerCropVolumeLogic::Apply(vtkMRMLCropVolumeParametersNode* pnode)
{
  vtkMRMLScene *scene = this->GetMRMLScene();

  vtkMRMLVolumeNode *inputVolume = 
    vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetInputVolumeNodeID()));
  vtkMRMLAnnotationROINode *inputROI = 
    vtkMRMLAnnotationROINode::SafeDownCast(scene->GetNodeByID(pnode->GetROINodeID()));

  if(!inputVolume || !inputROI)
    {
    std::cerr << "Failed to look up input volume and/or ROI!" << std::endl;
    return -1;
    }

  vtkMRMLScalarVolumeNode *refVolume;
  vtkMRMLVolumeNode *outputVolume = NULL;
  vtkMatrix4x4 *inputRASToIJK = vtkMatrix4x4::New();
  vtkMatrix4x4 *inputIJKToRAS = vtkMatrix4x4::New();
  vtkMatrix4x4 *outputRASToIJK = vtkMatrix4x4::New();
  vtkMatrix4x4 *outputIJKToRAS = vtkMatrix4x4::New();
  vtkMRMLLinearTransformNode *movingVolumeTransform = NULL, *roiTransform = NULL;

  // make sure inputs are initialized
  if(!inputVolume || !inputROI )
    {
    std::cerr << "CropVolume: Inputs are not initialized" << std::endl;
    return -1;
    }

  // check the output volume type
  vtkMRMLDiffusionTensorVolumeNode *dtvnode= vtkMRMLDiffusionTensorVolumeNode::SafeDownCast(inputVolume);
  vtkMRMLDiffusionWeightedVolumeNode *dwvnode= vtkMRMLDiffusionWeightedVolumeNode::SafeDownCast(inputVolume);
  vtkMRMLVectorVolumeNode *vvnode= vtkMRMLVectorVolumeNode::SafeDownCast(inputVolume);
  vtkMRMLScalarVolumeNode *svnode = vtkMRMLScalarVolumeNode::SafeDownCast(inputVolume);

  if(!this->Internal->VolumesLogic){
      std::cerr << "CropVolume: ERROR: failed to get hold of Volumes logic" << std::endl;
      return -2;
  }
  if(dtvnode){
    std::cerr << "CropVolume: ERROR: Diffusion tensor volumes are not supported by this module!" << std::endl;
    return -2;
  }
 
  std::ostringstream outSS;
  double outputSpacing[3], spacingScaleConst = pnode->GetSpacingScalingConst();
  outSS << inputVolume->GetName() << "-subvolume-scale_" << spacingScaleConst;

 if(dwvnode){
    outputVolume = vtkMRMLVolumeNode::SafeDownCast(this->Internal->VolumesLogic->CloneVolume(this->GetMRMLScene(), inputVolume, outSS.str().c_str()));
  }
  if(vvnode){
    outputVolume = vtkMRMLVolumeNode::SafeDownCast(this->Internal->VolumesLogic->CloneVolume(this->GetMRMLScene(), inputVolume, outSS.str().c_str()));
  }
  if(svnode){
    outputVolume = vtkMRMLVolumeNode::SafeDownCast(this->Internal->VolumesLogic->CloneVolume(this->GetMRMLScene(), inputVolume, outSS.str().c_str()));
  }
  refVolume = this->Internal->VolumesLogic->CreateLabelVolume(this->GetMRMLScene(), inputVolume, "CropVolume_ref_volume");
  refVolume->HideFromEditorsOn();

  //vtkMatrix4x4 *volumeXform = vtkMatrix4x4::New();
  //vtkMatrix4x4 *roiXform = vtkMatrix4x4::New();
  //vtkMatrix4x4 *T = vtkMatrix4x4::New();

  refVolume->GetRASToIJKMatrix(inputRASToIJK);
  refVolume->GetIJKToRASMatrix(inputIJKToRAS);
  outputRASToIJK->Identity();
  outputIJKToRAS->Identity();

  //T->Identity();
  //roiXform->Identity();
  //volumeXform->Identity();

  // prepare the resampling reference volume
  double roiRadius[3], roiXYZ[3];
  inputROI->GetRadiusXYZ(roiRadius);
  inputROI->GetXYZ(roiXYZ);
  std::cerr << "ROI radius: " << roiRadius[0] << "," << roiRadius[1] << "," << roiRadius[2] << std::endl;
  std::cerr << "ROI center: " << roiXYZ[0] << "," << roiXYZ[1] << "," << roiXYZ[2] << std::endl;

  double* inputSpacing = inputVolume->GetSpacing();
  double minSpacing = inputSpacing[0];
  if (minSpacing > inputSpacing[1])
    {
    minSpacing = inputSpacing[1];
    }
  if (minSpacing > inputSpacing[2])
    {
    minSpacing = inputSpacing[2];
    }

  if(pnode->GetIsotropicResampling()){
      outputSpacing[0] = minSpacing*spacingScaleConst;
      outputSpacing[1] = minSpacing*spacingScaleConst;
      outputSpacing[2] = minSpacing*spacingScaleConst;
  } else {
      outputSpacing[0] = inputSpacing[0]*spacingScaleConst;
      outputSpacing[1] = inputSpacing[1]*spacingScaleConst;
      outputSpacing[2] = inputSpacing[2]*spacingScaleConst;
  }

  int outputExtent[3];

  outputExtent[0] = roiRadius[0]/outputSpacing[0]*2.;
  outputExtent[1] = roiRadius[1]/outputSpacing[1]*2.;
  outputExtent[2] = roiRadius[2]/outputSpacing[2]*2.;

  outputIJKToRAS->SetElement(0,0,outputSpacing[0]);
  outputIJKToRAS->SetElement(1,1,outputSpacing[1]);
  outputIJKToRAS->SetElement(2,2,outputSpacing[2]);

  outputIJKToRAS->SetElement(0,3,roiXYZ[0]-roiRadius[0]+outputSpacing[0]*.5);
  outputIJKToRAS->SetElement(1,3,roiXYZ[1]-roiRadius[1]+outputSpacing[1]*.5);
  outputIJKToRAS->SetElement(2,3,roiXYZ[2]-roiRadius[2]+outputSpacing[2]*.5);

  // account for the ROI parent transform, if present
  roiTransform = vtkMRMLLinearTransformNode::SafeDownCast(inputROI->GetParentTransformNode());
  if(roiTransform){
    vtkMatrix4x4 *roiMatrix = vtkMatrix4x4::New();
    roiTransform->GetMatrixTransformToWorld(roiMatrix);
    outputIJKToRAS->Multiply4x4(roiMatrix, outputIJKToRAS, outputIJKToRAS);
  }

  outputRASToIJK->DeepCopy(outputIJKToRAS);
  outputRASToIJK->Invert();

  vtkImageData* outputImageData = vtkImageData::New();
  outputImageData->SetDimensions(outputExtent[0], outputExtent[1], outputExtent[2]);
  outputImageData->AllocateScalars();

  refVolume->SetAndObserveImageData(outputImageData);
  outputImageData->Delete();

  refVolume->SetIJKToRASMatrix(outputIJKToRAS);
  refVolume->SetRASToIJKMatrix(outputRASToIJK);

  inputRASToIJK->Delete();
  inputIJKToRAS->Delete();
  outputRASToIJK->Delete();
  outputIJKToRAS->Delete();

  if(this->Internal->ResampleLogic == 0)
    {
    std::cerr << "CropVolume: ERROR: resample logic is not set!";
    return -3;
    }

  vtkSmartPointer<vtkMRMLCommandLineModuleNode> cmdNode =
    this->Internal->ResampleLogic->CreateNodeInScene();
  assert(cmdNode.GetPointer() != 0);

  cmdNode->SetParameterAsString("inputVolume", inputVolume->GetID());
  cmdNode->SetParameterAsString("referenceVolume",refVolume->GetID());
  cmdNode->SetParameterAsString("outputVolume",outputVolume->GetID());

  movingVolumeTransform = vtkMRMLLinearTransformNode::SafeDownCast(inputVolume->GetParentTransformNode());

  if(movingVolumeTransform != NULL)
    cmdNode->SetParameterAsString("transformationFile",movingVolumeTransform->GetID());

  std::string interp = "linear";
  switch(pnode->GetInterpolationMode()){
    case 1: interp = "nn"; break;
    case 2: interp = "linear"; break;
    case 3: interp = "ws"; break;
    case 4: interp = "bs"; break;
  }

  cmdNode->SetParameterAsString("interpolationType", interp.c_str());
  this->Internal->ResampleLogic->ApplyAndWait(cmdNode);

  this->GetMRMLScene()->RemoveNode(refVolume);
  this->GetMRMLScene()->RemoveNode(cmdNode);

  outputVolume->SetAndObserveTransformNodeID(NULL);
  pnode->SetOutputVolumeNodeID(outputVolume->GetID());

  return 0;
}


//----------------------------------------------------------------------------
void vtkSlicerCropVolumeLogic::CropVoxelBased(vtkMRMLAnnotationROINode* roi, vtkMRMLVolumeNode* inputVolume, vtkMRMLVolumeNode* outputVolume)
{
  if(!roi || !inputVolume || !outputVolume)
    return;

  vtkNew<vtkImageData> imageDataWorkingCopy;
  imageDataWorkingCopy->DeepCopy(inputVolume->GetImageData());

  vtkNew<vtkMatrix4x4> inputRASToIJK;
  inputVolume->GetRASToIJKMatrix(inputRASToIJK.GetPointer());

  vtkNew<vtkMatrix4x4> inputIJKToRAS;
  inputVolume->GetIJKToRASMatrix(inputIJKToRAS.GetPointer());

  double roiXYZ[3];
  double roiRadius[3];

  roi->GetXYZ(roiXYZ);
  roi->GetRadiusXYZ(roiRadius);

  std::cout << "ROI XYZ: X="<<roiXYZ[0]<<" Y="<<roiXYZ[1]<<" Z="<<roiXYZ[2] << std::endl;
  std::cout << "RADIUS XYZ: X="<<roiRadius[0]<<" Y="<<roiRadius[1]<<" Z="<<roiRadius[2] << std::endl;

  double minXYZRAS[] = {roiXYZ[0]-roiRadius[0], roiXYZ[1]-roiRadius[1],roiXYZ[2]-roiRadius[2],1.};
  double maxXYZRAS[] = {roiXYZ[0]+roiRadius[0], roiXYZ[1]+roiRadius[1],roiXYZ[2]+roiRadius[2],1.};

  std::cout << "MinCorner in RAS: X="<<minXYZRAS[0]<<" Y="<<minXYZRAS[1]<<" Z="<<minXYZRAS[2] << std::endl;
  std::cout << "MaxCorner in RAS: X="<<maxXYZRAS[0]<<" Y="<<maxXYZRAS[1]<<" Z="<<maxXYZRAS[2] << std::endl;

  double minXYZIJK[4], maxXYZIJK[4];

  //transform to ijk
  inputRASToIJK->MultiplyPoint(minXYZRAS, minXYZIJK);
  inputRASToIJK->MultiplyPoint(maxXYZRAS, maxXYZIJK);

  std::cout << "MinCorner in IJK: X="<<minXYZIJK[0]<<" Y="<<minXYZIJK[1]<<" Z="<<minXYZIJK[2] << std::endl;
  std::cout << "MaxCorner in IJK: X="<<maxXYZIJK[0]<<" Y="<<maxXYZIJK[1]<<" Z="<<maxXYZIJK[2] << std::endl;

  int outputWholeExtent[6] = {std::min(minXYZIJK[0]+0.5,maxXYZIJK[0]+0.5),std::max(minXYZIJK[0]+0.5,maxXYZIJK[0]+0.5),std::min(minXYZIJK[1]+0.5,maxXYZIJK[1]+0.5),std::max(minXYZIJK[1]+0.5,maxXYZIJK[1]+0.5),std::min(minXYZIJK[2]+0.5,maxXYZIJK[2]+0.5),std::max(minXYZIJK[2]+0.5,maxXYZIJK[2]+0.5)};

  std::cout << "Whole Extent X values: Xmin="<<outputWholeExtent[0]<<" Xmax="<<outputWholeExtent[1]<< std::endl;
  std::cout << "Whole Extent Y values: Ymin="<<outputWholeExtent[2]<<" Ymax="<<outputWholeExtent[3]<< std::endl;
  std::cout << "Whole Extent Z values: Zmin="<<outputWholeExtent[4]<<" Zmax="<<outputWholeExtent[5]<< std::endl;


  double ijkNewOrigin[] = {outputWholeExtent[0],outputWholeExtent[2],outputWholeExtent[4],1.0};
  double  rasNewOrigin[4];
  inputIJKToRAS->MultiplyPoint(ijkNewOrigin,rasNewOrigin);

  vtkNew<vtkImageClip> imageClip;
  imageClip->SetInput(imageDataWorkingCopy.GetPointer());
  imageClip->SetOutputWholeExtent(outputWholeExtent);
  imageClip->SetClipData(true);

  imageClip->Update();

  vtkNew<vtkMatrix4x4> outputIJKToRAS;
  outputIJKToRAS->DeepCopy(inputIJKToRAS.GetPointer());

  outputIJKToRAS->SetElement(0,3,rasNewOrigin[0]);
  outputIJKToRAS->SetElement(1,3,rasNewOrigin[1]);
  outputIJKToRAS->SetElement(2,3,rasNewOrigin[2]);

  outputVolume->SetOrigin(rasNewOrigin[0],rasNewOrigin[1],rasNewOrigin[2]);


 // const char* scanOrder = vtkMRMLVolumeNode::ComputeScanOrderFromIJKToRAS(inputIJKToRAS.GetPointer());
 // bool computeIJKToRAS = vtkMRMLVolumeNode::ComputeIJKToRASFromScanOrder(scanOrder,inputVolume->GetSpacing(),imageClip->GetOutput()->GetDimensions(),false,outputIJKToRAS.GetPointer());


  vtkNew<vtkMatrix4x4> outputRASToIJK;
  outputRASToIJK->DeepCopy(outputIJKToRAS.GetPointer());
  outputRASToIJK->Invert();

  vtkNew<vtkImageData> outputImageData;
  outputImageData->DeepCopy(imageClip->GetOutput());

  int extent[6];
  imageClip->GetOutput()->GetExtent(extent);

  outputImageData->SetExtent(0, extent[1]-extent[0]-1, 0, extent[3]-extent[2]-1, 0, extent[5]-extent[4]-1);

  outputVolume->SetAndObserveImageData(outputImageData.GetPointer());

  outputVolume->SetIJKToRASMatrix(outputIJKToRAS.GetPointer());
  outputVolume->SetRASToIJKMatrix(outputRASToIJK.GetPointer());

  outputVolume->Modified();

}

//----------------------------------------------------------------------------
void vtkSlicerCropVolumeLogic::RegisterNodes()
{
  if(!this->GetMRMLScene())
    {
    return;
    }
  vtkMRMLCropVolumeParametersNode* pNode = vtkMRMLCropVolumeParametersNode::New();
  this->GetMRMLScene()->RegisterNodeClass(pNode);
  pNode->Delete();
}
