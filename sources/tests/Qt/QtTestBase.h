/**
 * \file QtTestBase.h
 * 2012-01-17 LB Initial implementation
 */

#ifndef QTTESTBASE_H
#define QTTESTBASE_H

#include <QObject>
#include <QTest>

class QString;

/// @brief This class should be used as the base class for Qt / Vtk tests.
/// It is based on QtTest framework. It allows to compare a string to a reference file.
class QtTestBase : public QObject
{
	Q_OBJECT 

public:
	/// @brief Constructor, if argc > 1 see compareToReference()
	QtTestBase(int argc);
	
	/// @brief Compares a string to the contents of a reference file.
	/// On differences it also writes an html of the diff to
	/// [build_dir]/tests/results/[refFile].html
	/// If the constructor is initialized with argc > 1 then reference file is
	/// updated with the contents of string.
	void compareToReference(QString string, QString refFile);
	
	/// @brief Returns the absolute path to the testdata input directory.
	QString getTestdataInputDir();

private:
	void writeReferenceFile(QString string, QString filePath);
	bool _writeRef;
};

#endif // QTTESTBASE_H