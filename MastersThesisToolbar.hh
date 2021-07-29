#pragma once

#if QT_VERSION >= 0x050000
#include <QtWidgets>
#else
#include <QtGui>
#endif

#include "ui_MastersThesisToolbarBase.h"

class MastersThesisToolbar : public QWidget, public Ui::MastersThesisToolbarBase
{
    Q_OBJECT

public:
    MastersThesisToolbar(QWidget * parent = 0);
};
