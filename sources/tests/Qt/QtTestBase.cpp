/**
 * \file QtTestBase.cpp
 * 2012-01-17 LB Implementation of QtTestBase class
 */

// ** INCLUDES **
#include "QtTestBase.h"

#include "Configure.h"

#include <QFile>
#include <QString>
#include <QTest>
#include <QTextStream>
#include <QRegExp>

#include "gdiff.h" // This should be the last include.

#include <iostream>

QtTestBase::QtTestBase(int argc) : _writeRef(false)
{
	if(argc > 1)
		_writeRef = true; 
}

void QtTestBase::compareToReference(QString string, QString refFile)
{
	QString filePath(TESTDATAPATH);
	filePath += "/ref/";
	filePath += refFile;
	
	QFile qFile(filePath);
	if(!qFile.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		if(_writeRef)
		{
			writeReferenceFile(string, filePath);
			return;
		}
		else
			QFAIL(QString("Reference file %1 could not be read.").arg(refFile).toAscii());
	}
		
	QString fileContent = qFile.readAll();
		
	// File compare
	diff_match_patch gdiff;
	QList<Diff> diffs = gdiff.diff_main(string, fileContent);
	
	bool equal = true;
	foreach(Diff diff, diffs)
	{
		if(diff.operation != EQUAL)
		{
			equal = false;
			break;
		}
	}

	if(!equal)
	{
		// Update reference file
		if(_writeRef)
		{
			writeReferenceFile(string, filePath);
			return;
		}
		
		// Html output
		QString htmlOutput = gdiff.diff_prettyHtml(diffs);
		QString htmlOutputFilename(BUILDPATH);
		htmlOutputFilename += "/tests/results/";
		htmlOutputFilename += refFile;
		htmlOutputFilename += ".html";
		QFile htmlFile(htmlOutputFilename);
		if (!htmlFile.open(QIODevice::WriteOnly | QIODevice::Text))
			QFAIL("Html output file could not be written.");
		QTextStream htmlStream(&htmlFile);
		htmlStream << htmlOutput;

#ifdef JENKINS_URL
		// On jenkins construct a link to the diff html file
		QString jenkinsURL(JENKINS_URL);
		QString jobName(JENKINS_JOB_NAME);
		QString relativePath;

		QRegExp rx("(.*)(/.*build[^/]*.*)");
		int pos = rx.indexIn(htmlOutputFilename);
		if (pos > -1)
			relativePath = rx.cap(2);

		// This assumes that the build dir is inside the sources dir.
		htmlOutputFilename = QString("%1job/%2/lastCompletedBuild/artifact/sources%3")
			.arg(jenkinsURL).arg(jobName).arg(relativePath);
#endif
		QFAIL(QString("File compare failed. See %1 for the differences.")
			.arg(htmlOutputFilename).toAscii());
	}
}

QString QtTestBase::getTestdataInputDir()
{
	QString dir(TESTDATAPATH);
	dir += "/input/";
	return dir;
}

void QtTestBase::writeReferenceFile(QString string, QString filePath)
{
	QFile file(filePath);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
		QWARN("Reference file could not be updated.");

	QTextStream out(&file);
	out << string;
	QWARN(QString("Reference file %1 updated.").arg(filePath).toAscii());
}