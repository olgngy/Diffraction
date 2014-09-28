#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <queue>
#include <QWidget>
#include <QMessageBox>
#include <QImage>
#include <QPainter>
#include <QLabel>
#include <QPixmap>
#include <QLineEdit>
#include <QPushButton>
#include "laserenv.h"
#include "utils.h"
#include "geomutils.h"
using namespace std;

// For x in [0, n)
void LaserEnv::fftrec(vector<cc>& a, bool inv) {
   int n = a.size();
   if (n == 1) return;
   vector<cc> a0(n/2), a1(n/2);
   for (int i = 0, j = 0; i < n; i+=2, j++)
      a0[j] = a[i], a1[j] = a[i+1];
   fftrec(a0, inv);
   fftrec(a1, inv);
   double ang = 2*Pi/n * (inv ? -1 : 1);
   cc w(1), wn(cos(ang), sin(ang));
   for (int i = 0; i < n/2; i++) {
      a[i] = a0[i] + w*a1[i];
      a[i+n/2] = a0[i] - w*a1[i];
      if (inv) a[i] /= 2, a[i+n/2] /= 2;
      w *= wn;
   }
}

// For x in [-n/2, n/2), y in [-n/2, n/2)
table LaserEnv::fourierTransform2D(const table& t, bool invert) {
   table ret(n, vector<cc>(n));
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         ret[i][j] = cc(0);
         for (int i2 = 0; i2 < n; i2++)
            for (int j2 = 0; j2 < n; j2++) {
               double ang1 = 2*Pi*(i-n/2)*(i2-n/2)/n * (invert ? -1 : 1);
               double ang2 = 2*Pi*(j-n/2)*(j2-n/2)/n * (invert ? -1 : 1);
               ret[i][j] += t[i2][j2] * cc(cos(ang1), sin(ang1)) * cc(cos(ang2), sin(ang2));
            }
         ret[i][j] /= n;
      }
   return ret;
}

// For x in [-n/2, n/2), y in [-n/2, n/2)
table LaserEnv::fastFourierTransform2D(const table& t, bool invert) {
   table T = t;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         if (invert) T[i][j] *= n;
         else T[i][j] /= n; // Adjust
   if (invert) // Swap diagonal quadrants
      for (int i = 0; i < n/2; i++)
         for (int j = 0; j < n/2; j++) {
            swap(T[i][j], T[i+n/2][j+n/2]);
            swap(T[i][j+n/2], T[i+n/2][j]);
         }
   for (int i = 0; i < n; i++)
      fftrec(T[i], invert);
   for (int i = 0; i < n; i++) // Transpose
      for (int j = i+1; j < n; j++)
         swap(T[i][j], T[j][i]);
   for (int i = 0; i < n; i++)
      fftrec(T[i], invert);
   for (int i = 0; i < n; i++) // Transpose
      for (int j = i+1; j < n; j++)
         swap(T[i][j], T[j][i]);
   if (!invert) // Swap diagonal quadrants
      for (int i = 0; i < n/2; i++)
         for (int j = 0; j < n/2; j++) {
            swap(T[i][j], T[i+n/2][j+n/2]);
            swap(T[i][j+n/2], T[i+n/2][j]);
         }
   return T;
}

void LaserEnv::print(const table& t) {
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++)
         printf("(%+.3lf %+.3lf)\t\t", t[i][j].real(), t[i][j].imag());
   puts("");
}

table LaserEnv::generateBaseRay() {
   table ret(n, vector<cc>(n));
   //double U0 = 1;
   double R0 = w*w/2;
   double w1 = w;
   double R1 = 0;
   double z1;
   double Phi;
   double maxi = 0;
//   if (z == 0) {
//       w1 = w;
//       R1 = 0;
//       Phi = 0;
//   }
//   else {
//       z1 = z;//Transformation();
//       w1 = z1*w*sqrt(1 + z1*z1/R0/R0);
//       R1 = z1 + R0*R0/z1;
//   }
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         double val = pow(exp(-((i-n/2.)*(i-n/2.) + (j-n/2.)*(j-n/2.))/(2.*w1*w1)), 2) * 1000;
         maxi = max(maxi, val);
      }
   vector<double> phiRow(n);
   for (int i = 0; i < n; i++) {
      //generatePhiForBaseRayRow(phiRow, i, w1, R1);
      for (int j = 0; j < n; j++) {
         Phi = phiRow[j];
         //double H = sqrt(R1*R1 - n*n/4);
         //double X = -n/2.0 + j;
         //double Y = sqrt(R1*R1 - X*X) - H;
         //Phi = Y/(2*Pi);
         //Phi = 200*Pi - j*200*Pi/1023;
         //Phi = 4*100*Pi*j*(1023 - j)/1023/1023;
         double val = pow(exp(-((i-n/2.)*(i-n/2.) + (j-n/2.)*(j-n/2.))/(2.*w1*w1)), 2) * 1000;
         ret[i][j] = cc(val/maxi*cos(Phi), val/maxi*sin(Phi));
      }
   }
   return ret;
}

double LaserEnv::zTransformation() {
    double w0 = w*0.001/40;
    double R0 = w0*w0/2;
    double wz = w0*sqrt(1 + z*z/(R0*R0));
    std::cout << w*wz/(w0*2) << endl;
    return w*wz/(w0*2);

}

void LaserEnv::generatePhiForBaseRayRow(vector<double>& phi, int row, double rCircle, double rCurve) {
   int l = -1, r = -1;
   for (int i = 0; i < n; i++) {
      double distFromCenter = (i-n/2.0)*(i-n/2.0) + (row-n/2.0)*(row-n/2.0);
      if (distFromCenter <= rCircle*rCircle + 1e-8) {
         // Inside the circle.
         if (l == -1) l = i;
         r = i;
      }
   }
   phi.assign(n, 0);
   if (l == -1) return;
   int width = r-l+1;
   double h = sqrt(rCurve*rCurve - width*width/4);
//   int DBG = 511;
//   if (row == DBG) {
//      cout << "==================================" << endl;
//      cout << "rCircle = " << rCircle << endl;
//      cout << "rCurve = " << rCurve << endl;
//      cout << "l = " << l << endl;
//      cout << "r = " << r << endl;
//      cout << "width = " << width << endl;
//   }
   for (int i = 0; i < width; i++) {
      double X = -width/2.0 + i;
      double Y = sqrt(rCurve*rCurve - X*X) - h;
      phi[l+i] = Y/(2*Pi);
   }
//   if (row == DBG) {
//     for (int i = 0; i < n; i++) cout << phi[i] << " ";
//     cout << endl;
//   }
}

void LaserEnv::processCircleCenter(QImage& img, int I, int J, vector<int>& circ) {
   int R = int(sz/3. + rand()*1./RAND_MAX*sz*2./3 + 0.5);
   int delta = sz/3;
   for (int i = 0; i < circ.size(); i+=3) {
      int I2 = circ[i];
      int J2 = circ[i+1];
      int R2 = circ[i+2];
      if ((I2-I)*(I2-I) + (J2-J)*(J2-J) <= (R+R2+delta)*(R+R2+delta))
         return;
   }
   circ.push_back(I);
   circ.push_back(J);
   circ.push_back(R);
}

QImage LaserEnv::drawCircularLattice(int canvasWidth) {
   QImage ret = QImage(canvasWidth, n, QImage::Format_Mono);
   for (int i = 0; i < n; i++)
      for (int j = 0; j < canvasWidth; j++) {
         ret.setPixel(j, i, 0);
         if (j < 0 || i < 0) {
             printf("Problem 1: I = %d   J = %d \n", i, j);
         }
      }
   int R = sz;
   int delta = 4*sz;
   int I = 2*R + delta;
   int k = 0;
   int x = 0;
   QPainter painter(&ret);
   painter.setBrush(QBrush(Qt::white));
   for (int i = 0; i < n; i+=I) {
       int ini = 0;
       if (k%2 == 1)
           ini = - R - delta;
       for (int j = ini; j < canvasWidth; j+=I) {
         painter.drawEllipse(j, i, 2*R, 2*R);
         x = rand()%(4*sz);
         j += x;
         //cout << delta + x << endl;
       }
       k++;
   }

//   vector<int> p1(n), p2(canvasWidth), circ;
//   for (int i = 0; i < n; i++)
//      p1[i] = i;
//   for (int i = 0; i < canvasWidth; i++)
//      p2[i] = i;
//   //random_shuffle(p1.begin(), p1.end());
//   //random_shuffle(p2.begin(), p2.end());
//   for (int i = 0; i < n; i++)
//      for (int j = 0; j < canvasWidth; j++)
//         processCircleCenter(ret, p1[i], p2[j], circ);
//   QPainter painter(&ret);
//   painter.setBrush(QBrush(Qt::white));
//   for (int i = 0; i < (int)circ.size(); i+=3) {
//      int I = circ[i], J = circ[i+1], R = circ[i+2];
//      painter.drawEllipse(J-R, I-R, 2*R, 2*R);
//   }
   return ret;
}

table LaserEnv::multiplyTables(const table& t1, const table& t2) {
   table ret(n, vector<cc>(n));
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         ret[i][j] = t1[i][j] * t2[i][j];
   return ret;
}

double LaserEnv::calculateRayPhi(double CIRCLE_CENTER_X,
                                 double CIRCLE_CENTER_Y,
                                 double CIRCLE_RADIUS,
                                 double RAY_Y) {
    const double N1 = 1.0;
    const double N2 = 1.5;
    pt ray1 = pt(0, RAY_Y), ray2 = pt(1, RAY_Y);
    circle c = circle(CIRCLE_CENTER_X, CIRCLE_CENTER_Y, CIRCLE_RADIUS);
    line l1 = line(ray1, ray2);
    vector<pt> intersectionPoints = line_circle_intersect(l1, c);
    //if (intersectionPoints.size() < 2)
    //    throw 1;
    pt pInter = intersectionPoints.front();
    pt pCenter = pt(CIRCLE_CENTER_X, CIRCLE_CENTER_Y);
    pt dirNormal = (pCenter-pInter).normalize();
    pt dirRay = (ray2-ray1).normalize();
    double theta1 = rotationAngle(dirNormal*-1, dirRay*-1);
    double theta2 = asin(sin(theta1)*N1/N2) + Pi;
    pt dirRefracted = (rotate(pInter-dirNormal, theta2, pInter)-pInter).normalize();
    line l2 = line(pInter, pInter+dirRefracted);
    vector<pt> intersectionPoints2 = line_circle_intersect(l2, c);
    pt pInter2 = intersectionPoints2.back();
    pt dirNormal2 = (pInter2-pCenter).normalize();
    double angleTouch = rotationAngle(dirNormal2*-1, dirRefracted*-1);
    pt dirRefracted2 = (rotate(pInter2-dirNormal2, -angleTouch, pInter2)-pInter2).normalize();
    dirRefracted2 = dirRefracted2*-1;
    line l3 = line(pInter2, pInter2+dirRefracted2);
    vector<pt> intersectionPoints3 = line_circle_intersect(l3, c);
    pt pInter3 = intersectionPoints3.front();
    pt dirNormal3 = (pInter3-pCenter).normalize();
    pt dirRay3 = (pInter3-pInter2).normalize();
    double theta3 = rotationAngle(dirNormal3, dirRay3);
    double theta4 = asin(sin(theta3)*N2/N1) + Pi;
    pt dirRefracted4 = (rotate(pInter3-dirNormal3, theta4, pInter3)-pInter3).normalize();
    line l4 = line(pInter3, pInter3+dirRefracted4);
    pt pFinish;
    line_line_intersect(l4, line(pt(0, 0), pt(0, 1)), pFinish);
    double L = (pInter2-pInter).len() + (pInter3-pInter2).len();
    const double Lambda = 0.6;
    double Phi = 2*L*N2 / Lambda;
    // Adjust to [0..2Pi)
    Phi = Phi - (int)(Phi/2)*2;
    Phi *= Pi;
    //cerr << "L = " << L << endl;
    return Phi;
}

table LaserEnv::read(const QImage& img) {
   table ret(n, vector<cc>(n));
   const double Pi = acos(0.0)*2;
   vector< vector<bool> > used(n, vector<bool>(n));
   const int di[] = {-1, 0, 1, 0};
   const int dj[] = {0, 1, 0, -1};
   queue<int> q;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         QRgb rgb = img.pixel(j, i);
         bool black = !qRed(rgb) && !qGreen(rgb) && !qBlue(rgb);
         if (black || used[i][j]) continue;
         vector<int> connected;
         q.push(i);
         q.push(j);
         connected.push_back(i);
         connected.push_back(j);
         used[i][j] = true;
         int mini = 1e9, maxi = -1e9;
         int minj = 1e9, maxj = -1e9;
         while (!q.empty()) {
            int vi = q.front(); q.pop();
            int vj = q.front(); q.pop();
            for (int k = 0; k < 4; k++) {
               int ni = vi+di[k];
               int nj = vj+dj[k];
               if (ni >= 0 && ni < n && nj >= 0 && nj < n && !used[ni][nj]) {
                  QRgb rgb2 = img.pixel(nj, ni);
                  bool black2 = !qRed(rgb2) && !qGreen(rgb2) && !qBlue(rgb2);
                  if (black2) continue;
                  used[ni][nj] = true;
                  q.push(ni);
                  q.push(nj);
                  connected.push_back(ni);
                  connected.push_back(nj);
                  mini = min(mini, ni);
                  maxi = max(maxi, ni);
                  minj = min(minj, nj);
                  maxj = max(maxj, nj);
               }
            }
         }
         double r = (1+rand())*1./(RAND_MAX+1.);
         //double phi = (1+rand())*1./(RAND_MAX+1.);
         //double standardNormalDist = cos(2*Pi*phi) * sqrt(-2*log(r));
         //double normalVal = deviation * standardNormalDist;
         double cx = 0.5*(mini+maxi);
         double cy = 0.5*(minj+maxj);
         double cr = max(maxi-mini, maxj-minj)+1;
         for (int k = 0; k < (int)connected.size(); k+=2) {
            // Point inside the circle
            int vi = connected[k];
            int vj = connected[k+1];
            double distToCenterOfCircle = (pt(cx, cy)-pt(vi, vj)).len();
            double ry = cy + distToCenterOfCircle;
            //cerr << "cx = " << cx << endl;
            //cerr << "cy = " << cy << endl;
            //cerr << "cr = " << cr << endl;
            double phiVal = calculateRayPhi(cx, cy, cr, ry);
            //cerr << "phiVal = " << phiVal << endl;
            //cerr << "normalVal = " << normalVal << endl;
            ret[vi][vj] = cc(cos(phiVal), sin(phiVal));
         }
      }
   return ret;
}

QImage LaserEnv::write(const table& t) {
   QImage ret = QImage(n, n, QImage::Format_Indexed8);
   QVector<QRgb> colors;
   for (int i = 0; i < 256; i++)
      colors.push_back(qRgb(i, i, i));
   ret.setColorTable(colors);
   double mini = 1e10, maxi = 0;
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
         double a = abs(t[i][j]);
         double v = log(a*a + 10*deviation);
         mini = min(mini, v);
         maxi = max(maxi, v);
      }
   //printf("N = %d\n", n);
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++) {
         double a = abs(t[i][j]);
         double v = log(a*a + 10*deviation);
         int col = ((v-mini)/(maxi-mini))*255;
         if (col < 0)
             col = 0;
         ret.setPixel(j, i, col);
         if (j < 0 || i < 0 || col < 0) {
             printf("Problem 2: I = %d   J = %d  Col = %d \n", i, j, col);
         }
         //printf("%d ", col);
      }
   return ret;
}

QLabel* LaserEnv::createLabel(const string& text) {
    QLabel* label = new QLabel(text.c_str(), this);
    label->setFont(QFont("Arial", 8));
    return label;
}

void LaserEnv::createGUIControls() {
    setWindowTitle("Author: Olga Nad");
    this->setFixedSize(1100, 365);
    const int d = 15;
    int n1 = 256;
    lblRay = new QLabel(this);
    lblRay->move(d, d);
    lblLattice = new QLabel(this);
    lblLattice->move(d+n1+d, d);
    lblDiffraction = new QLabel(this);
    lblDiffraction->move(d+(n1+d)*2, d);
    lblFFT = new QLabel(this);
    lblFFT->move(d+(n1+d)*3, d);
    lblLatticeText = createLabel("Lattice");
    lblLatticeText->move(d+n1/3+n1+d, d+n1+3);
    lblLatticeText->show();
    lblRayText = createLabel("Laser beam");
    lblRayText->move(d+n1/3, d+n1+3);
    lblRayText->show();
    lblDiffractionText = createLabel("Interference");
    lblDiffractionText->move(d+n1/3+(n1+d)*2, d+n1+3);
    lblDiffractionText->show();
    lblFFTText = createLabel("Diffraction");
    lblFFTText->move(d+n1/3+(n1+d)*3, d+n1+3);
    lblFFTText->show();
    lblSize = createLabel("Size of lattice (px):");
    lblSize->move(d, n1+d*2+3);
    lblSize->show();
    lblWidth = createLabel("Width of beam (px):");
    lblWidth->move(d, d/2+n1+d*3+1);
    lblWidth->show();
    lblDeviation = createLabel("Deviation:");
    lblDeviation->move(d+300, n1+d*2+3);
    lblDeviation->show();
    lblIterations = createLabel("Iterations:");
    lblIterations->move(d+300, d/2+n1+d*3+1);
    lblIterations->show();
    lblZ = createLabel("Z(cm):");
    lblZ->move(d+180, n1+d*2+3);
    lblZ->show();
    lblFrameX = createLabel("Frame X-coord (px):");
    lblFrameX->move(d+540, n1+d*2+3);
    lblFrameX->show();
    lblFrameY = createLabel("Frame Y-coord (px):");
    lblFrameY->move(d+540, d/2+n1+d*3+1);
    lblFrameY->show();
    lblFrameWidth = createLabel("Frame width (px):");
    lblFrameWidth->move(d+730, n1+d*2+3);
    lblFrameWidth->show();
    lblFrameHeight = createLabel("Frame height (px):");
    lblFrameHeight->move(d+730, d/2+n1+d*3+1);
    lblFrameHeight->show();
    lblMotion = createLabel("Velocity (px):");
    lblMotion->move(d+910, n1+d*2+3);
    lblMotion->show();
    lblIterSpectrum = createLabel("Spectrum Iterations:");
    lblIterSpectrum->move(d+910, d/2+n1+d*3+3);
    lblIterSpectrum->show();
    editSize = new QLineEdit(this);
    editSize->move(d+130, n1+d*2+1);
    editSize->resize(40, 20);
    editSize->show();
    editWidth = new QLineEdit(this);
    editWidth->move(d+130, d/2+n1+d*3+1);
    editWidth->resize(40, 20);
    editWidth->show();
    editDeviation = new QLineEdit(this);
    editDeviation->move(d+370, n1+d*2+1);
    editDeviation->resize(40, 20);
    editDeviation->show();
    editIterations = new QLineEdit(this);
    editIterations->move(d+370, d/2+n1+d*3+1);
    editIterations->resize(40, 20);
    editIterations->show();
    editZ = new QLineEdit(this);
    editZ->move(d+250, n1+d*2+1);
    editZ->resize(40, 20);
    editZ->show();
    editFrameX = new QLineEdit(this);
    editFrameX->move(d+680, n1+d*2+1);
    editFrameX->resize(40, 20);
    editFrameX->show();
    editFrameY = new QLineEdit(this);
    editFrameY->move(d+680, d/2+n1+d*3+1);
    editFrameY->resize(40, 20);
    editFrameY->show();
    editFrameHeight = new QLineEdit(this);
    editFrameHeight->move(d+860, d/2+n1+d*3+1);
    editFrameHeight->resize(40, 20);
    editFrameHeight->show();
    editFrameWidth = new QLineEdit(this);
    editFrameWidth->move(d+860, n1+d*2+1);
    editFrameWidth->resize(40, 20);
    editFrameWidth->show();
    editMotion = new QLineEdit(this);
    editMotion->move(d+1000, n1+d*2+1);
    editMotion->resize(40, 20);
    editMotion->show();
    editIterSpectrum = new QLineEdit(this);
    editIterSpectrum->move(d+1000, d/2+n1+d*3+1);
    editIterSpectrum->resize(40, 20);
    editIterSpectrum->show();

    generateButton = new QPushButton(this);
    generateButton->setText("Generate");
    generateButton->resize(110, 25);
    generateButton->move(d+180, d+n1+d*4-2);
    generateButton->show();

    iterateButton = new QPushButton(this);
    iterateButton->setText("Print average");
    iterateButton->resize(115, 25);
    iterateButton->move(d+420, d+n1+d*4-2);
    iterateButton->show();

    motionButton = new QPushButton(this);
    motionButton->setText("Move");
    motionButton->resize(135, 25);
    motionButton->move(d+910, d+n1+d*4-2);
    motionButton->show();
    connect(generateButton, SIGNAL(clicked()), this, SLOT(generateClicked()));
    connect(iterateButton, SIGNAL(clicked()), this, SLOT(iterateClicked()));
    connect(motionButton, SIGNAL(clicked()), this, SLOT(motionClicked()));
}

void LaserEnv::setupGUIControls() {
   // Laser beam
   table tabRay = generateBaseRay();
   QImage imgRay = write(tabRay);
   int lblW = 256;
   int lblH = 256;
   lblRay->setPixmap(QPixmap::fromImage(imgRay).scaled(lblW, lblH, Qt::KeepAspectRatio));
   lblRay->show();
   imgRay.save("Ray.bmp", "bmp");
   // Circular lattice
   QImage imgLattice = drawCircularLattice(n);
   lblLattice->setPixmap(QPixmap::fromImage(imgLattice).scaled(lblW, lblH, Qt::KeepAspectRatio));
   lblLattice->show();
   imgLattice.save("Lattice.bmp", "bmp");
   // Diffraction
   table tabLattice = read(imgLattice);
   table tabDiffraction = multiplyTables(tabRay, tabLattice);
   QImage imgDiffraction = write(tabDiffraction);
   lblDiffraction->setPixmap(QPixmap::fromImage(imgDiffraction).scaled(lblW, lblH, Qt::KeepAspectRatio));
   lblDiffraction->show();
   imgDiffraction.save("Diffraction.bmp", "bmp");
   // Diffraction's FFT
   table t4 = fastFourierTransform2D(tabDiffraction, false);
   QImage imgFFT = write(t4);
   lblFFT->setPixmap(QPixmap::fromImage(imgFFT).scaled(lblW, lblH, Qt::KeepAspectRatio));
   lblFFT->show();
   imgFFT.save("DiffractionFFT.bmp", "bmp");
   // Write text output
   /*freopen("diffraction.txt", "wt", stdout);
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++)
         printf("%.2lf %.2lf\t", t4[i][j].real(), t4[i][j].imag());
   freopen("diffraction_abs.txt", "wt", stdout);
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++)
         printf("%.2lf, ", abs(t4[i][j]));
   */
   // Control panel
   editSize->setText(intToString(sz).c_str());
   editWidth->setText(intToString(w).c_str());
   editDeviation->setText(doubleToString(deviation).c_str());
   editIterations->setText(intToString(averageIterations).c_str());
   editZ->setText(doubleToString(z).c_str());
   editFrameX->setText(intToString(frameX).c_str());
   editFrameY->setText(intToString(frameY).c_str());
   editFrameWidth->setText(intToString(frameWidth).c_str());
   editFrameHeight->setText(intToString(frameHeight).c_str());
   editMotion->setText(intToString(motionSpeed).c_str());
   editIterSpectrum->setText(intToString(iterSpectrum).c_str());
}

void LaserEnv::generateClicked() {
   sz = atoi(editSize->text().toStdString().c_str());
   w = atoi(editWidth->text().toStdString().c_str());
   z = atoi(editZ->text().toStdString().c_str());
   averageIterations = atoi(editIterations->text().toStdString().c_str());
   motionSpeed = atoi(editMotion->text().toStdString().c_str());
   frameX = atoi(editFrameX->text().toStdString().c_str());
   frameY = atoi(editFrameY->text().toStdString().c_str());
   frameHeight = atoi(editFrameHeight->text().toStdString().c_str());
   frameWidth = atoi(editFrameWidth->text().toStdString().c_str());
   iterSpectrum = atoi(editIterSpectrum->text().toStdString().c_str());
   sscanf(editDeviation->text().toStdString().c_str(), "%lf", &deviation);
   setupGUIControls();
}

void LaserEnv::iterateClicked() {
    int lblW = 256;
    int lblH = 256;
   generateClicked();
   table tabRay = generateBaseRay();
   QImage imgRay = write(tabRay);
   QImage imgLattice = drawCircularLattice(n);
   table tDiffFFTAverage(n, vector<cc>(n));
   for (int i = 0; i < averageIterations; i++) {
      table tabLattice = read(imgLattice);
      table tabDiffraction = multiplyTables(tabRay, tabLattice);
      table tFFT = fastFourierTransform2D(tabDiffraction, false);
      for (int j = 0; j < n; j++)
         for (int k = 0; k < n; k++)
            tDiffFFTAverage[j][k] += cc(abs(tFFT[j][k]), 0);
      //QImage imgFFT = write(tFFT);
      //char fileName[64];
      //sprintf(fileName, "DiffractionDFT-%d.png", i+1);
      //imgFFT.save(fileName, "png");
   }
   for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
         tDiffFFTAverage[j][k] /= averageIterations;
   QImage imgDiffFFTAverage = write(tDiffFFTAverage);
   imgDiffFFTAverage.save("DiffractionDFT-AVERAGE.png", "png");
   QLabel* imageLabel = new QLabel();
   imageLabel->setPixmap(QPixmap::fromImage(imgDiffFFTAverage));
   imageLabel->setFixedSize(imgDiffFFTAverage.width(), imgDiffFFTAverage.height());
   imageLabel->show();
   // Write text output
   freopen("diffraction_average.txt", "wt", stdout);
   for (int i = 0; i < n; i++, putchar('\n'))
      for (int j = 0; j < n; j++)
         printf("%.2lf %.2lf\t", tDiffFFTAverage[i][j].real(), tDiffFFTAverage[i][j].imag());
   freopen("diffraction_average_abs.txt", "wt", stdout);
   for (int i = 0; i < n; i++) {
      printf("{");
      for (int j = 0; j < n; j++) {
         printf("%.2lf", abs(tDiffFFTAverage[i][j]));
         if (j != n-1) printf(", ");
      }
      printf("},\n");
   }
   //QMessageBox msgBox;
   //msgBox.setText("Resultant image generated!");
   //msgBox.exec();
}

//QImage LaserEnv::drawCircularLattice2(int canvasWidth, int i1) {
//    circleNumber = 0;
//   QImage ret = QImage(canvasWidth, n, QImage::Format_Mono);
//   for (int i = 0; i < n; i++)
//      for (int j = 0; j < canvasWidth; j++)
//         ret.setPixel(j, i, 0);
//   QPainter painter(&ret);
//   painter.setBrush(QBrush(Qt::white));
//      int I = n/2, J = i1, R = 20/2;
//      painter.drawEllipse(J-R, I-R, 2*R, 2*R);
//      /*if (I > frameX && I < frameX + frameWidth && J > frameY && J < frameY + frameHeight)
//          circleNumber++;
//          */
//   //cout << circleNumber << endl;
//   return ret;
//}

void LaserEnv::motionClicked() {
   generateClicked();
   table tabRay = generateBaseRay();
   QImage imgRay = write(tabRay);
   int len = n + n;
   QImage imgLattice;
   drawCircularLattice(len);
   char fileName[64];
   vector<double> meanValuesIter;
   for (int iter = 0; iter < iterSpectrum; iter++) {
        vector<cc> meanValues;
        imgLattice = drawCircularLattice(len);
        imgLattice.save("Lattice1.bmp", "bmp");
        //int i11 = 20;

        for (int shift = 0, it = 0; shift < len - n; shift += motionSpeed, it++) {
            QImage curLattice = QImage(n, n, QImage::Format_Mono);
            //QImage curLattice = drawCircularLattice2(1024, i11);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    QRgb rgb = imgLattice.pixel(j+shift, i);
                    bool black = !qRed(rgb) && !qGreen(rgb) && !qBlue(rgb);
                    curLattice.setPixel(j, i, !black);
                    if (j < 0 || i < 0) {
                        printf("Problem 3: I = %d   J = %d \n", i, j);
                    }
                }
            table tabLattice = read(curLattice);
            table tabDiffraction = multiplyTables(tabRay, tabLattice);
            table tabFFT = fastFourierTransform2D(tabDiffraction, false);
            QImage writeLattice = write(tabDiffraction);
            QImage writeFFT = write(tabFFT);
            const int i1 = frameY, j1 = frameX;
            const int i2 = frameY+frameHeight, j2 = frameX+frameWidth;
        //      for (int i = i1; i <= i2; i++) {
        //         writeFFT.setPixel(j1, i, 255);
        //         writeFFT.setPixel(j2, i, 255);
        //      }
        //      for (int j = j1; j <= j2; j++) {
        //         writeFFT.setPixel(j, i1, 255);
        //         writeFFT.setPixel(j, i2, 255);
        //      }
              cc curMean;
              int cnt = 0;
              for (int i = i1; i <= i2; i++)
                 for (int j = j1; j <= j2; j++)
                    curMean += tabFFT[i][j], cnt++;
              curMean /= cnt;
              meanValues.push_back(curMean);
        //      char fileName[64];
        //      sprintf(fileName, "CurLattice-%d.png", it+1);
        //      curLattice.save(fileName, "png");
        //      sprintf(fileName, "Motion-%d.png", it+1);
        //      writeLattice.save(fileName, "png");
        //      sprintf(fileName, "Motion-Diffraction-%d.png", it+1);
        //      writeFFT.save(fileName, "png");
        //      i11++;
        }

    //   printf("Initial values: ");
    //   double MV = 0;
    //   for (int i = 0; i < meanValues.size(); i++) {
    //      double meanAbs = abs(meanValues[i]);
    //      MV += meanAbs;
    //      printf("%.12lf\t", meanAbs);
    //   }
    //   fflush(stdout);
    //   QImage prespectrumImage = writeSpectrumImage(meanValues);
    //   prespectrumImage.save("!Motion-PRESPECTRUM.png", "png");
    //   QLabel* imageLabel1 = new QLabel();
    //   imageLabel1->setPixmap(QPixmap::fromImage(prespectrumImage));
    //   imageLabel1->setFixedSize(prespectrumImage.width(), prespectrumImage.height());
    //   imageLabel1->show();
       /*
       puts("");
       MV /= meanValues.size();
       double DV = 0;
       for (int i = 0; i < meanValues.size(); i++) {
          double meanAbs = abs(meanValues[i]);
          DV += (meanAbs-MV)*(meanAbs-MV);
       }
       DV = sqrt(DV);
       printf("Mean = %.12lf\n", MV);
       printf("Standard deviation = %.12lf\n", DV);
       printf("Standard deviation / Mean = %.12lf\n", DV/MV);
       printf("Mean / Standard deviation = %.12lf\n", MV/DV);
       */
        sprintf(fileName, "Motion-%d.txt", iter);
        freopen(fileName, "wt", stdout);
        for (int i = 0; i < meanValues.size(); i++)
           printf("%.12lf\n", abs(meanValues[i]));
        puts("");
        fflush(stdout);

        fftrec(meanValues, false);
       for (int j = 0; j < meanValues.size()/2; j++) {
          swap(meanValues[j], meanValues[j+meanValues.size()/2]);

       }
       if (iter == 0 ) {
           for (int i = 0; i < meanValues.size(); i++)
               meanValuesIter.push_back(abs(meanValues[i]));
       }
       else {
           for (int i = 0; i < meanValues.size(); i++)
               meanValuesIter[i] += abs(meanValues[i]);
       }
       sprintf(fileName, "MotionDFT-%d.txt", iter);
       freopen(fileName, "wt", stdout);
       for (int i = 0; i < meanValues.size(); i++)
          printf("%.12lf\n", abs(meanValues[i]));
       puts("");
       fflush(stdout);
   }

   for (int i = 0; i < meanValuesIter.size(); i++) {
       meanValuesIter[i] /= iterSpectrum;
   }
   double maxValue = *max_element(meanValuesIter.begin(), meanValuesIter.end());
   for (int i = 0; i < meanValuesIter.size(); i++) {
       meanValuesIter[i] /= maxValue;
   }

   sprintf(fileName, "MotionDFTResult.txt");
   freopen(fileName, "wt", stdout);
   for (int i = 0; i < meanValuesIter.size(); i++)
      printf("%.12lf\n", meanValuesIter[i]);
   puts("");
   fflush(stdout);

//   QImage spectrumImage = writeSpectrumImage(meanValues);
//   spectrumImage.save("!Motion-SPECTRUM.png", "png");
   QLabel* imageLabel = new QLabel();
//   imageLabel->setPixmap(QPixmap::fromImage(spectrumImage));
//   imageLabel->setFixedSize(spectrumImage.width(), spectrumImage.height());
   imageLabel->show();
}

QImage LaserEnv::writeSpectrumImage(const vector<cc>& row) {
   int len = row.size();
   double mini = 1e10, maxi = -1e10;
   for (int i = 0; i < len; i++) {
      double val = abs(row[i]);
      mini = min(mini, val);
      maxi = max(maxi, val);
   }
   int step = 5;
   int height = 200;
   int width = (len+1)*step;
   QImage retImage = QImage(width, height, QImage::Format_Mono);
   for (int i = 0; i < width; i++)
      for (int j = 0; j < height; j++)
         retImage.setPixel(i, j, 0);
   for (int i = 0; i < width; i++)
      retImage.setPixel(i, height-10, 1);
   for (int i = 0; i < len; i++) {
      double val = (abs(row[i]) - mini) / (maxi - mini);
      int h = val * (height-10);
      for (int j = 0; j < h; j++)
         retImage.setPixel(step*(i+1), height-j-10, 1);
   }
   return retImage;
}

LaserEnv::LaserEnv(QWidget* obj)
   : QWidget(obj), n(1024)
{
   sz = 15;
   w = 15;
   z = 1;
   deviation = 1;
   averageIterations = 5;
   motionSpeed = 1;
   frameX = 511;
   frameY = 511;
   frameHeight = 512;
   frameWidth = 512;
   iterSpectrum = 1;
   createGUIControls();
   setupGUIControls();
}
