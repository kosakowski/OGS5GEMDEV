/**
 * \file TestHelperFunctions.cpp
 * 2012-01-20 LB Implementation of TestHelperFunctions
 */

// ** INCLUDES **
#include "TestHelperFunctions.h"

#include "Configure.h"

#include <QFile>
#include <QString>
#include <QTest>
#include <QTextStream>
#include <QRegExp>

#include "gtest.h"
#include "gdiff.h" // This should be the last include.

#include <iostream>

void compareToReference(QString string, QString refFile)
{
	QString filePath(TESTDATAPATH);
	filePath += "/ref/";
	filePath += refFile;
	
	QFile qFile(filePath);
	if(!qFile.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		if(g_writeRef)
		{
			writeReferenceFile(string, filePath);
			return;
		}
		else
		{
			std::cout << QString("Reference file %1 could not be read.").arg(refFile).toStdString() << std::endl;
			ASSERT_TRUE(false);
		}
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
		if(g_writeRef)
		{
			writeReferenceFile(string, filePath);
			return;
		}
		
		// Html output
		QString htmlOutput = gdiff.diff_prettyHtml(diffs);
		QString htmlOutputFilename(BUILDPATH);
		htmlOutputFilename += "/tests/results/";
		htmlOutputFilename += refFile;
		if(g_writeOutput)
			writeOutputFile(string, htmlOutputFilename);
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
		std::cout << QString("File compare failed. See %1 for the differences.")
			.arg(htmlOutputFilename).toStdString() << std::endl;
		ASSERT_TRUE(false);
	}
}

QString getTestdataInputDir()
{
	QString dir(TESTDATAPATH);
	dir += "/input/";
	return dir;
}

void writeReferenceFile(QString string, QString filePath)
{
	QFile file(filePath);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		std::cout << "Reference file could not be updated." << std::endl;
		return;
	}

	QTextStream out(&file);
	out << string;
	std::cout << QString("Reference file %1 updated.").arg(filePath).toStdString() << std::endl;
}

void writeOutputFile(QString string, QString filePath)
{
	QFile file(filePath);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		std::cout << "Output file could not be written." << std::endl;
		return;
	}

	QTextStream out(&file);
	out << string;
	std::cout << QString("Output file %1 updated.").arg(filePath).toStdString() << std::endl;
}
