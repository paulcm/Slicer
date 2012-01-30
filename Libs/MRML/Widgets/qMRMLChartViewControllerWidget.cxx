/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Qt includes
#include <QActionGroup>
#include <QDebug>
#include <QInputDialog>
#include <QLabel>
#include <QMenu>
#include <QHBoxLayout>

// CTK includes
#include <ctkLogger.h>
#include <ctkPopupWidget.h>

// qMRML includes
#include "qMRMLColors.h"
#include "qMRMLNodeFactory.h"
#include "qMRMLSceneViewMenu.h"
#include "qMRMLChartView.h"
#include "qMRMLChartViewControllerWidget_p.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLChartViewNode.h>
#include <vtkMRMLChartNode.h>
#include <vtkMRMLSceneViewNode.h>

// STD include
#include <string>

//--------------------------------------------------------------------------
static ctkLogger logger("org.slicer.libs.qmrmlwidgets.qMRMLChartViewControllerWidget");
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// qMRMLChartViewControllerWidgetPrivate methods

//---------------------------------------------------------------------------
qMRMLChartViewControllerWidgetPrivate::qMRMLChartViewControllerWidgetPrivate(
  qMRMLChartViewControllerWidget& object)
  : Superclass(object)
{
  this->ChartViewNode = 0;
  this->ChartView = 0;
}

//---------------------------------------------------------------------------
qMRMLChartViewControllerWidgetPrivate::~qMRMLChartViewControllerWidgetPrivate()
{
}

//---------------------------------------------------------------------------
void qMRMLChartViewControllerWidgetPrivate::setupPopupUi()
{
  Q_Q(qMRMLChartViewControllerWidget);

  this->Superclass::setupPopupUi();
  this->PopupWidget->setAlignment(Qt::AlignBottom | Qt::AlignLeft);
  this->Ui_qMRMLChartViewControllerWidget::setupUi(this->PopupWidget);
  
  // configure the Chart selector
  this->chartComboBox->addAttribute("vtkMRMLChartNode", "Chart", "1");

  // Connect Chart selector
  this->connect(this->chartComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
                SLOT(onChartNodeSelected(vtkMRMLNode*)));

  // configure the Array selector
  this->arrayComboBox->addAttribute("vtkMRMLDoubleArrayNode", "Array", "1");

  // Connect the actions
  QObject::connect(this->actionShow_Lines, SIGNAL(toggled(bool)),
                   q, SLOT(showLines(bool)));
  QObject::connect(this->actionShow_Markers, SIGNAL(toggled(bool)),
                   q, SLOT(showMarkers(bool)));
  QObject::connect(this->actionShow_Grid, SIGNAL(toggled(bool)),
                   q, SLOT(showGrid(bool)));
  QObject::connect(this->actionShow_Legend, SIGNAL(toggled(bool)),
                   q, SLOT(showLegend(bool)));
  
  // Connect the buttons
  this->showLinesToolButton->setDefaultAction(this->actionShow_Lines);
  this->showMarkersToolButton->setDefaultAction(this->actionShow_Markers);
  this->showGridToolButton->setDefaultAction(this->actionShow_Grid);
  this->showLegendToolButton->setDefaultAction(this->actionShow_Legend);

  // Connect the checkboxes
  QObject::connect(this->showTitleCheckBox, SIGNAL(toggled(bool)),
                   q, SLOT(showTitle(bool)));
  QObject::connect(this->showXAxisLabelCheckBox, SIGNAL(toggled(bool)),
                   q, SLOT(showXAxisLabel(bool)));
  QObject::connect(this->showYAxisLabelCheckBox, SIGNAL(toggled(bool)),
                   q, SLOT(showYAxisLabel(bool)));

  // Connect the text boxes
  QObject::connect(this->titleLineEdit, SIGNAL(textEdited(const QString&)),
                   q, SLOT(setTitle(const QString&)));
  QObject::connect(this->xAxisLabelLineEdit, SIGNAL(textEdited(const QString&)),
                   q, SLOT(setXAxisLabel(const QString&)));
  QObject::connect(this->yAxisLabelLineEdit, SIGNAL(textEdited(const QString&)),
                   q, SLOT(setYAxisLabel(const QString&)));
  
  // Connect the scene
  QObject::connect(q, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)),
                   this->chartComboBox, SLOT(setMRMLScene(vtkMRMLScene*)));
  
}

//---------------------------------------------------------------------------
void qMRMLChartViewControllerWidgetPrivate::init()
{
  this->Superclass::init();
  this->ViewLabel->setText(qMRMLChartViewControllerWidget::tr("1"));
  this->BarLayout->addStretch(1);
  //this->setColor(QColor("#6e4b7c"));
  this->setColor(QColor("#e1ba3c"));
}


// --------------------------------------------------------------------------
vtkMRMLChartNode* qMRMLChartViewControllerWidgetPrivate::chartNode()
{
  Q_Q(qMRMLChartViewControllerWidget);

  if (!this->ChartViewNode || !q->mrmlScene())
    {
    // qDebug() << "No ChartViewNode or no Scene";
    return 0;
    }

  // Get the current chart node
  vtkMRMLChartNode *chartNode 
    = vtkMRMLChartNode::SafeDownCast(q->mrmlScene()->GetNodeByID(this->ChartViewNode->GetChartNodeID()));

  return chartNode;
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidgetPrivate::onChartNodeSelected(vtkMRMLNode * node)
{
  Q_Q(qMRMLChartViewControllerWidget);

  if (!this->ChartViewNode)
    {
    return;
    }

  if (this->chartNode() == node)
    {
    return;
    }

  this->qvtkReconnect(this->chartNode(), node, vtkCommand::ModifiedEvent,
                      q, SLOT(updateWidgetFromMRML()));
  
  this->ChartViewNode->SetChartNodeID(node ? node->GetID() : 0);

  if (node)
    {
    q->updateWidgetFromMRML();
    }
}



// --------------------------------------------------------------------------
// qMRMLChartViewControllerWidget methods

// --------------------------------------------------------------------------
qMRMLChartViewControllerWidget::qMRMLChartViewControllerWidget(QWidget* parentWidget)
  : Superclass(new qMRMLChartViewControllerWidgetPrivate(*this), parentWidget)
{
  Q_D(qMRMLChartViewControllerWidget);
  d->init();
}

// --------------------------------------------------------------------------
qMRMLChartViewControllerWidget::~qMRMLChartViewControllerWidget()
{
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::setChartView(qMRMLChartView* view)
{
  Q_D(qMRMLChartViewControllerWidget);
  d->ChartView = view;
}

//---------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::setViewLabel(const QString& newViewLabel)
{
  Q_D(qMRMLChartViewControllerWidget);

  if (d->ChartViewNode)
    {
    logger.error("setViewLabel should be called before setViewNode !");
    return;
    }

  d->ChartViewLabel = newViewLabel;
  d->ViewLabel->setText(d->ChartViewLabel);
}

//---------------------------------------------------------------------------
CTK_GET_CPP(qMRMLChartViewControllerWidget, QString, viewLabel, ChartViewLabel);


// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::setMRMLChartViewNode(
    vtkMRMLChartViewNode * viewNode)
{
  Q_D(qMRMLChartViewControllerWidget);
  this->qvtkReconnect(d->ChartViewNode, viewNode, vtkCommand::ModifiedEvent,
                      this, SLOT(updateWidgetFromMRML()));
  d->ChartViewNode = viewNode;
  this->updateWidgetFromMRML();
}

// --------------------------------------------------------------------chart------
void qMRMLChartViewControllerWidget::updateWidgetFromMRML()
{
  Q_D(qMRMLChartViewControllerWidget);
  
  //qDebug() << "qMRMLChartViewControllerWidget::updateWidgetFromMRML()";

  if (!d->ChartViewNode || !this->mrmlScene())
    {
    return;
    }

  vtkMRMLChartNode *chartNode = d->chartNode();
  if (!chartNode)
    {
    return;
    }

  // ChartNode selector
  d->chartComboBox->setCurrentNode(chartNode->GetID());
  
  // Buttons
  const char *propertyValue;
  propertyValue = chartNode->GetProperty("default", "showLines");
  d->actionShow_Lines->setChecked(propertyValue && !strcmp("on", propertyValue));

  propertyValue = chartNode->GetProperty("default", "showMarkers");
  d->actionShow_Markers->setChecked(propertyValue && !strcmp("on", propertyValue));

  propertyValue = chartNode->GetProperty("default", "showGrid");
  d->actionShow_Grid->setChecked(propertyValue && !strcmp("on", propertyValue));

  propertyValue = chartNode->GetProperty("default", "showLegend");
  d->actionShow_Legend->setChecked(propertyValue && !strcmp("on", propertyValue));

  // Titles, axis labels (checkboxes AND text widgets)
  propertyValue = chartNode->GetProperty("default", "showTitle");
  d->showTitleCheckBox->setChecked(propertyValue && !strcmp("on", propertyValue));
  propertyValue = chartNode->GetProperty("default", "title");
  if (propertyValue)
    {
    d->titleLineEdit->setText(propertyValue);
    }
  else
    {
    d->titleLineEdit->clear();
    }

  propertyValue = chartNode->GetProperty("default", "showXAxisLabel");
  d->showXAxisLabelCheckBox->setChecked(propertyValue && !strcmp("on", propertyValue));
  propertyValue = chartNode->GetProperty("default", "xAxisLabel");
  if (propertyValue)
    {
    d->xAxisLabelLineEdit->setText(propertyValue);
      }
  else
    {
    d->xAxisLabelLineEdit->clear();
    }

  propertyValue = chartNode->GetProperty("default", "showYAxisLabel");
  d->showYAxisLabelCheckBox->setChecked(propertyValue && !strcmp("on", propertyValue));
  propertyValue = chartNode->GetProperty("default", "yAxisLabel");
  if (propertyValue)
    {
    d->yAxisLabelLineEdit->setText(propertyValue);
    }
  else
    {
    d->yAxisLabelLineEdit->clear();
    }

}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::setMRMLScene(vtkMRMLScene* newScene)
{
  Q_D(qMRMLChartViewControllerWidget);
  
  qDebug() << "Inside setMRMLScene()";

  if (this->mrmlScene() == newScene)
    {
    return;
    }

  d->qvtkReconnect(this->mrmlScene(), newScene, vtkMRMLScene::EndBatchProcessEvent,
                   this, SLOT(updateWidgetFromMRML()));

  // Disable the node selectors as they would fire signal currentIndexChanged(0)
  // meaning that there is no current node anymore. It's not true, it just means 
  // that the current node was not in the combo box list menu before
  bool chartBlockSignals = d->chartComboBox->blockSignals(true);
  bool arrayBlockSignals = d->arrayComboBox->blockSignals(true);

  this->Superclass::setMRMLScene(newScene);

  d->chartComboBox->blockSignals(chartBlockSignals);
  d->arrayComboBox->blockSignals(arrayBlockSignals);

  if (this->mrmlScene())
    {
    this->updateWidgetFromMRML();
    }
}
  
// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::showLines(bool show)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "showLines", show ? "on" : "off");

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "showLines");
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::showMarkers(bool show)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "showMarkers", show ? "on" : "off");

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "showMarkers");
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::showGrid(bool show)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "showGrid", show ? "on" : "off");

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "showGrid");
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::showLegend(bool show)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "showLegend", show ? "on" : "off");

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "showLegend");
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::showTitle(bool show)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "showTitle", show ? "on" : "off");

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "showTitle");
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::showXAxisLabel(bool show)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "showXAxisLabel", show ? "on" : "off");

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "showXAxisLabel");
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::showYAxisLabel(bool show)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "showYAxisLabel", show ? "on" : "off");

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "showYAxisLabel");
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::setTitle(const QString &str)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "title", str.toStdString().c_str());

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "title");
}

// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::setXAxisLabel(const QString &str)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "xAxisLabel", str.toStdString().c_str());

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "xAxisLabel");
}


// --------------------------------------------------------------------------
void qMRMLChartViewControllerWidget::setYAxisLabel(const QString &str)
{
  Q_D(qMRMLChartViewControllerWidget);

  vtkMRMLChartNode *chartNode = d->chartNode();

  if (!chartNode)
    {
    return;
    }

  // Set the parameter
  chartNode->SetProperty("default", "yAxisLabel", str.toStdString().c_str());

  //qDebug() << "Regetting property: " << chartNode->GetProperty("default", "yAxisLabel");
}
