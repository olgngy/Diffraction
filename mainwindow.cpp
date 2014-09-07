#include "mainwindow.h"
#include "utils.h"
#include "geomutils.h"

#include <iostream>

#include <QPainter>
#include <QLabel>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent) {

    WINDOW_WIDTH = 1024;
    WINDOW_HEIGHT = 690;
    WINDOW_IDENT = 50;
    CIRCLE_CENTER_X = 500;
    CIRCLE_CENTER_Y = 300;
    CIRCLE_RADIUS = 200;
    RAY_Y = 200;
    N1 = 1.0;
    N2 = 1.5;

    setFixedSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    setWindowTitle("Olga Nad Inc.™ All rights reserved.");

    installUserInterface();
}

MainWindow::~MainWindow() {
}

void MainWindow::installUserInterface() {
    const int EH = 20;
    const int EW = 45;

    editCircleCenterX = new QLineEdit(this);
    editCircleCenterX->setText(intToString(CIRCLE_CENTER_X).c_str());
    editCircleCenterX->resize(EW, EH);
    editCircleCenterX->show();

    editCircleCenterY = new QLineEdit(this);
    editCircleCenterY->setText(intToString(CIRCLE_CENTER_Y).c_str());
    editCircleCenterY->resize(EW, EH);
    editCircleCenterY->show();

    editCircleRadius = new QLineEdit(this);
    editCircleRadius->setText(intToString(CIRCLE_RADIUS).c_str());
    editCircleRadius->resize(EW, EH);
    editCircleRadius->show();

    editRayY = new QLineEdit(this);
    editRayY->setText(intToString(RAY_Y).c_str());
    editRayY->resize(EW, EH);
    editRayY->show();

    editN1 = new QLineEdit(this);
    editN1->setText(doubleToString(N1).c_str());
    editN1->resize(EW, EH);
    editN1->show();

    editN2 = new QLineEdit(this);
    editN2->setText(doubleToString(N2).c_str());
    editN2->resize(EW, EH);
    editN2->show();

    updateButton = new QPushButton(this);
    updateButton->setText("Update");
    updateButton->resize(280, 60);
    updateButton->show();

    connect(updateButton, SIGNAL(clicked()), this, SLOT(updateButtonClicked()));
}

void MainWindow::paintEvent(QPaintEvent* paintEvent) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    drawPlot(&painter);
    drawRay(&painter);
    drawUserInterface(&painter);
}

void MainWindow::drawPlot(QPainter* painter) {
    const int W = WINDOW_WIDTH;
    const int H = WINDOW_HEIGHT;
    const int ID = WINDOW_IDENT;

    // Fill white window
    painter->fillRect(ID, ID, W-2*ID, H-3*ID, Qt::white);

    // Draw X axis
    painter->drawLine(W-ID-10, ID-10, W-ID, ID-1); // Arrows
    painter->drawLine(W-ID-10, ID+10, W-ID, ID+1); // Arrows
    painter->drawLine(ID-20, ID, W-ID, ID); // Axis
    for (int i = ID; i < H-2*ID; i += 50) {
        painter->drawLine(ID-3, i, ID+3, i); // Dashes
        if (i != ID)
            painter->drawText(ID-35, i-10, 50, 20, 0, intToString(i-ID).c_str());
    }

    // Draw Y axis
    painter->drawLine(ID-10, H-2*ID-10, ID-1, H-2*ID); // Arrows
    painter->drawLine(ID+10, H-2*ID-10, ID+1, H-2*ID); // Arrows
    painter->drawLine(ID, ID-20, ID, H-2*ID); // Axis
    for (int i = ID; i < W-3*ID; i += 50) {
        painter->drawLine(i, ID-3, i, ID+3); // Dashes
        if (i != ID)
            painter->drawText(i-15, ID-25, 50, 20, 0, intToString(i-ID).c_str());
    }

    // Draw zero label
    painter->drawText(ID-15, ID-20, 20, 20, 0, "0");

    // Draw circle
    QPen oldPen = painter->pen();
    painter->setPen(QPen(Qt::black, 2));
    painter->setBrush(QBrush(QColor(98, 245, 237, 0x80)));
    painter->drawEllipse(ID+CIRCLE_CENTER_X-CIRCLE_RADIUS, ID+CIRCLE_CENTER_Y-CIRCLE_RADIUS,
                         2*CIRCLE_RADIUS, 2*CIRCLE_RADIUS);
    painter->setPen(oldPen);
}

void MainWindow::drawRay(QPainter* painter) {
    const int W = WINDOW_WIDTH;
    const int H = WINDOW_HEIGHT;
    const int ID = WINDOW_IDENT;
    const double RY = RAY_Y;
    pt ray1 = pt(0, RY), ray2 = pt(W-2*ID, RY);

    circle c = circle(CIRCLE_CENTER_X, CIRCLE_CENTER_Y, CIRCLE_RADIUS);
    line l1 = line(ray1, ray2);
    vector<pt> intersectionPoints = line_circle_intersect(l1, c);

    // Initial ray
    if (intersectionPoints.size() < 2) {
        painter->drawLine(ID, ID+RY, W-ID, ID+RY);
        L = Phi = -1;
    }
    else {
        QPen oldPen = painter->pen();
        pt pInter = intersectionPoints.front();
        pt pCenter = pt(CIRCLE_CENTER_X, CIRCLE_CENTER_Y);
        painter->setPen(QPen(Qt::red, 3));
        // Draw red line (1st segment)
        painter->drawLine(ID, ID+RY, ID+pInter.x, ID+pInter.y);
        line lineNormal = line(pInter, pCenter);
        pt dirInterface = pt(lineNormal.a, lineNormal.b).normalize();
        pt pLeft = pInter - dirInterface*CIRCLE_RADIUS;
        pt pRight = pInter + dirInterface*CIRCLE_RADIUS;
        // Draw dashed interface line
        painter->setPen(Qt::DashLine);
        painter->drawLine(ID+pLeft.x, ID+pLeft.y, ID+pRight.x, ID+pRight.y);
        pt dirNormal = (pCenter-pInter).normalize();
        pt pTop = pInter - dirNormal*CIRCLE_RADIUS;
        pt pBottom = pInter + dirNormal*CIRCLE_RADIUS;
        // Draw dashed normal line
        painter->setPen(Qt::DashLine);
        painter->drawLine(ID+pTop.x, ID+pTop.y, ID+pBottom.x, ID+pBottom.y);
        pt dirRay = (ray2-ray1).normalize();
        double theta1 = rotationAngle(dirNormal*-1, dirRay*-1);
        double theta2 = asin(sin(theta1)*N1/N2) + Pi;
        pt dirRefracted = (rotate(pInter-dirNormal, theta2, pInter)-pInter).normalize();
        line l2 = line(pInter, pInter+dirRefracted);
        vector<pt> intersectionPoints2 = line_circle_intersect(l2, c);
        pt pInter2 = intersectionPoints2.back();
        // Draw red line (2nd segment)
        painter->setPen(QPen(Qt::red, 3));
        painter->drawLine(ID+pInter.x, ID+pInter.y, ID+pInter2.x, ID+pInter2.y);
        pt dirNormal2 = (pInter2-pCenter).normalize();
        pt pTop2 = pInter2 - dirNormal2*CIRCLE_RADIUS;
        pt pBottom2 = pInter2 + dirNormal2*CIRCLE_RADIUS;
        // Draw dashed normal line (2nd segment)
        painter->setPen(Qt::DashLine);
        painter->drawLine(ID+pTop2.x, ID+pTop2.y, ID+pBottom2.x, ID+pBottom2.y);
        line lineNormal2 = line(pInter2, pCenter);
        pt dirInterface2 = pt(lineNormal2.a, lineNormal2.b).normalize();
        pt pLeft2 = pInter2 - dirInterface2*CIRCLE_RADIUS;
        pt pRight2 = pInter2 + dirInterface2*CIRCLE_RADIUS;
        // Draw dashed interface line (2nd segment)
        painter->setPen(Qt::DashLine);
        painter->drawLine(ID+pLeft2.x, ID+pLeft2.y, ID+pRight2.x, ID+pRight2.y);
        double angleTouch = rotationAngle(dirNormal2*-1, dirRefracted*-1);
        pt dirRefracted2 = (rotate(pInter2-dirNormal2, -angleTouch, pInter2)-pInter2).normalize();
        dirRefracted2 = dirRefracted2*-1;
        line l3 = line(pInter2, pInter2+dirRefracted2);
        vector<pt> intersectionPoints3 = line_circle_intersect(l3, c);
        pt pInter3 = intersectionPoints3.front();
        // Draw red line (2nd segment)
        painter->setPen(QPen(Qt::red, 3));
        painter->drawLine(ID+pInter2.x, ID+pInter2.y, ID+pInter3.x, ID+pInter3.y);
        pt dirNormal3 = (pInter3-pCenter).normalize();
        // Draw dashed normal line (2nd segment)
        pt pTop3 = pInter3 - dirNormal3*CIRCLE_RADIUS;
        pt pBottom3 = pInter3 + dirNormal3*CIRCLE_RADIUS;
        painter->setPen(Qt::DashLine);
        painter->drawLine(ID+pTop3.x, ID+pTop3.y, ID+pBottom3.x, ID+pBottom3.y);
        line lineNormal3 = line(pCenter, pInter3);
        pt dirInterface3 = pt(lineNormal3.a, lineNormal3.b).normalize();
        pt pLeft3 = pInter3 - dirInterface3*CIRCLE_RADIUS;
        pt pRight3 = pInter3 + dirInterface3*CIRCLE_RADIUS;
        // Draw dashed interface line (3nd segment)
        painter->setPen(Qt::DashLine);
        painter->drawLine(ID+pLeft3.x, ID+pLeft3.y, ID+pRight3.x, ID+pRight3.y);
        pt dirRay3 = (pInter3-pInter2).normalize();
        double theta3 = rotationAngle(dirNormal3, dirRay3);
        double theta4 = asin(sin(theta3)*N2/N1) + Pi;
        pt dirRefracted4 = (rotate(pInter3-dirNormal3, theta4, pInter3)-pInter3).normalize();
        line l4 = line(pInter3, pInter3+dirRefracted4);
        pt pFinish;
        line_line_intersect(l4, line(pt(0, 0), pt(0, 1)), pFinish);
        // Draw red line (4th segment)
        painter->setPen(QPen(Qt::red, 3));
        painter->drawLine(ID+pInter3.x, ID+pInter3.y, ID+pFinish.x, ID+pFinish.y);
        painter->setPen(oldPen);
        L = (pInter2-pInter).len() + (pInter3-pInter2).len();
        const double Lambda = 0.6;
        Phi = 2*Pi*L*N2 / Lambda;
        Phi = Phi - (int)(Phi/2)*2;
        Phi *= Pi;
    }
}

void MainWindow::drawUserInterface(QPainter* painter) {
    const int W = WINDOW_WIDTH;
    const int H = WINDOW_HEIGHT;
    const int ID = WINDOW_IDENT;
    const double RY = RAY_Y;
    const int Y1 = H-1.3*ID;
    const int Y2 = H-0.7*ID;

    painter->drawText(ID, Y1, "Circle center X:");
    painter->drawText(ID, Y2, "Circle center Y:");
    painter->drawText(ID+180, Y1, "Circle radius:");
    painter->drawText(ID+224, Y2, "Ray Y:");
    painter->drawText(ID+363, Y1, "Coeff. N1:");
    painter->drawText(ID+363, Y2, "Coeff. N2:");
    painter->drawText(ID+525, Y1, ("L = " + doubleToString(L)).c_str());
    painter->drawText(ID+525, Y2, ("Phi = " + doubleToString(Phi)).c_str());

    editCircleCenterX->move(ID+107, Y1-14);
    editCircleCenterY->move(ID+107, Y2-14);
    editCircleRadius->move(ID+272, Y1-14);
    editRayY->move(ID+272, Y2-14);
    editN1->move(ID+437, Y1-14);
    editN2->move(ID+437, Y2-14);

    updateButton->move(ID+650, Y1-18);
}

void MainWindow::updateButtonClicked() {
    sscanf(editCircleCenterX->text().toStdString().c_str(), "%lf", &CIRCLE_CENTER_X);
    sscanf(editCircleCenterY->text().toStdString().c_str(), "%lf", &CIRCLE_CENTER_Y);
    sscanf(editCircleRadius->text().toStdString().c_str(), "%lf", &CIRCLE_RADIUS);
    sscanf(editRayY->text().toStdString().c_str(), "%lf", &RAY_Y);
    sscanf(editN1->text().toStdString().c_str(), "%lf", &N1);
    sscanf(editN2->text().toStdString().c_str(), "%lf", &N2);
    repaint();
}
