/**
 * \file TestHelperFunctions.h
 * 2012-01-20 LB Initial implementation
 */

#ifndef TESTHELPERFUNCTIONS_H
#define TESTHELPERFUNCTIONS_H

class QString;

extern bool g_writeRef;
extern bool g_writeOutput;

void compareToReference(QString string, QString refFile);

QString getTestdataInputDir();

void writeReferenceFile(QString string, QString filePath);

void writeOutputFile(QString string, QString filePath);

#endif // TESTHELPERFUNCTIONS_H