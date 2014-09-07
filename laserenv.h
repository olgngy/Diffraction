#ifndef LASERENV_H
#define LASERENV_H

#include "geomutils.h"

#include <complex>
#include <vector>
#include <QWidget>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
using namespace std;

typedef complex<double> cc;
typedef vector< vector<cc> > table;

class LaserEnv : public QWidget
{
   Q_OBJECT
   const int n;
   int circleNumber;
   int sz;
   int w;
   int z;
   int ray1X;
   int ray1Y;
   int ray2X;
   int ray2Y;
   double deviation;
   int averageIterations;
   int motionSpeed;
   int frameX;
   int frameY;
   int frameWidth;
   int frameHeight;
   int iterSpectrum;
   QLabel* lblRay;
   QLabel* lblRayText;
   QLabel* lblLattice;
   QLabel* lblLatticeText;
   QLabel* lblDiffraction;
   QLabel* lblDiffractionText;
   QLabel* lblFFT;
   QLabel* lblFFTText;
   QLabel* lblSize;
   QLabel* lblWidth;
   QLabel* lblZ;
   QLabel* lblDeviation;
   QLabel* lblIterations;
   QLabel* lblMotion;
   QLabel* lblFrameX;
   QLabel* lblFrameY;
   QLabel* lblFrameHeight;
   QLabel* lblFrameWidth;
   QLabel* lblIterSpectrum;
   QLineEdit* editSize;
   QLineEdit* editWidth;
   QLineEdit* editZ;
   QLineEdit* editDeviation;
   QLineEdit* editIterations;
   QLineEdit* editMotion;
   QLineEdit* editFrameX;
   QLineEdit* editFrameY;
   QLineEdit* editFrameHeight;
   QLineEdit* editFrameWidth;
   QLineEdit* editIterSpectrum;
   QPushButton* iterateButton;
   QPushButton* generateButton;
   QPushButton* motionButton;
private:
   void  fftrec(vector<cc>& a, bool inv);
   table fourierTransform2D(const table& t, bool invert);
   table fastFourierTransform2D(const table& t, bool invert);
   void  print(const table& t);
   table generateBaseRay();
   double zTransformation();
   void generatePhiForBaseRayRow(vector<double>& phi, int row, double rCircle, double rCurve);
   QImage drawCircularLattice(int canvasWidth);
   QImage drawCircularLattice2(int canvasWidth, int i1);
   void processCircleCenter(QImage& img, int I, int J, vector<int>& circ);
   table multiplyTables(const table& t1, const table& t2);
   table read(const QImage& img);
   QImage write(const table& t);
   QImage writeSpectrumImage(const vector<cc>& row);
   void setupGUIControls();
   void createGUIControls();
   QLabel* createLabel(const string& text);
   double calculateRayPhi(double, double, double, double);
public:
   LaserEnv(QWidget* = 0);
private slots:
   void iterateClicked();
   void generateClicked();
   void motionClicked();
};

#endif // LASERENV_H
