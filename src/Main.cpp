/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "Parameters.h"
#include "BgzfFileType.h"
#include "VcfFileReader.h"
#include "VcfFileWriter.h"
#include "StringBasics.h"

void vcfVersion()
{
    std::cerr << "Version: " << VERSION
              << "; Built: " << DATE << " by "<< USER << std::endl << std::endl;
}

void description()
{
    std::cerr << " vcfRefGen - Clean/Reduce VCF files removing the info field, saving only the GT genotype field, and removing any records where any kept sample is not phased or is missing a genotype" << std::endl;
}

void usage()
{
    vcfVersion();
    description();
    std::cerr << "\t./vcfRefGen --in <input VCF File> --out <output VCF File> [--uncompress] [--sample <sampleFile.txt>] [--minAC <minAlleleCount>] [--filterList <chr start end> [--params]"<< std::endl;
    std::cerr << "\tRequired Parameters:\n"
              << "\t\t--in      : VCF file to read\n"
              << "\t\t--out     : VCF file to write\n"
              << "\tOptional Parameters:\n"
              << "\t\t--uncompress   : write an uncompressed VCF output file\n"
              << "\t\t--sampleSubset : file with samples IDs to keep.\n"
              << "\t\t--minAC        : min minor allele count to keep\n"
              << "\t\t--filterList   : regions to include, format: chr start end\n"
              << "\t\t--params       : print the parameter settings\n"
              << std::endl;
}


int main(int argc, char ** argv)
{
    String refFile = "";
    String inputVcf = "";
    String outputVcf = "";
    String sampleSubset = "";
    String filterList = "";
    int minAC = -1;
    bool uncompress = false;
    bool params = false;
    
    // Read in the parameters.    
    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inputVcf)
        LONG_STRINGPARAMETER("out", &outputVcf)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_PARAMETER("uncompress", &uncompress)
        LONG_STRINGPARAMETER("sampleSubset", &sampleSubset)
        LONG_INTPARAMETER("minAC", &minAC)
        LONG_STRINGPARAMETER("filterList", &filterList)
        LONG_PARAMETER("params", &params)
        END_LONG_PARAMETERS();

    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));
    
    inputParameters.Read(argc, &(argv[0]));
    
    // Check that all files were specified.
    if(inputVcf == "")
    {
        usage();
        inputParameters.Status();
        std::cerr << "Missing \"--in\", a required parameter.\n\n";
        return(-1);
    }
    if(outputVcf == "")
    {
        usage();
        inputParameters.Status();
        std::cerr << "Missing \"--out\", a required parameter.\n\n";
        return(-1);
    }

    if(params)
    {
        inputParameters.Status();
    }

    VcfFileReader inFile;
    VcfFileWriter outFile;
    VcfHeader header;
    VcfRecord record;

    // Open the file.
    if(sampleSubset.IsEmpty())
    {
        inFile.open(inputVcf, header);        
    }
    else
    {
        inFile.open(inputVcf, header, sampleSubset, NULL, NULL);
    }
    
    if(uncompress)
    {
        outFile.open(outputVcf, header, InputFile::DEFAULT);
    }
    else
    {
        outFile.open(outputVcf, header);
    }

    // Add the discard rule for minor allele count.
    if(minAC >= 0)
    {
        inFile.addDiscardMinMinorAlleleCount(minAC, NULL);
    }

    // Do not keep any records with a missing GT or with a non-phased sample.
    inFile.addDiscardRules(VcfFileReader::DISCARD_NON_PHASED | 
                           VcfFileReader::DISCARD_MISSING_GT);

    int numReadRecords = 0;
    int numWrittenRecords = 0;
    int returnVal = 0;

    // Set to only store/write the GT field.
    VcfRecordGenotype::addStoreField("GT");
    while(inFile.readRecord(record))
    {
        ++numReadRecords;
        
        // Clear the INFO field.
        record.getInfo().clear();
        // Write the record.
        if(!outFile.writeRecord(record))
        {
            // Write error.
            std::cerr << "Failed writing a vcf record.\n";
            returnVal = -1;
        }
        ++numWrittenRecords;
    }
 
    std::cerr << "NumReadRecords: " << inFile.getNumRecords()
              << "; NumRecordsDiscarded: "
              << inFile.getNumRecords() - numReadRecords
              << "; NumWrittenRecords: " << numWrittenRecords << "\n";

    inFile.close();   
    outFile.close();   

    return(returnVal);
}
