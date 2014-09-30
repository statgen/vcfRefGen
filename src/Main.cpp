/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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
#include "IntervalTree.h"

void vcfVersion()
{
    std::cerr << "Version: " << VERSION
              << "; Built: " << DATE << " by "<< USER << std::endl << std::endl;
}

void description()
{
    std::cerr << " vcfRefGen - Generate VCF reference panel for imputation" << std::endl;
    std::cerr << "             Clean/Reduce VCF files removing the info field, saving only the GT genotype field,\n"
              << "             and removing any records where any kept sample is not phased or is missing a genotype\n\n" << std::endl;
}

void usage()
{
    vcfVersion();
    description();
    std::cerr << "./vcfRefGen --in <input VCF File> --out <output VCF File> [--uncompress] [--sample <sampleFile.txt>] [--minAC <minAlleleCount>] [--filterList <filterFile> [--keepGT <field1,field2>] [--params]"<< std::endl;
    std::cerr << "\tRequired Parameters:\n"
              << "\t\t--in      : VCF file to read\n"
              << "\t\t--out     : VCF file to write\n"
              << "\tOptional Parameters:\n"
              << "\t\t--allfields    : keep info & all genotype fields and\n"
              << "\t\t                 do not filter out non-phased or missing genotype\n"
              << "\t\t                 genotype records.\n"
              << "\t\t--uncompress   : write an uncompressed VCF output file\n"
              << "\t\t--sampleSubset : file with samples IDs to keep.\n"
              << "\t\t--minAC        : min minor allele count to keep\n"
              << "\t\t--filterList   : filename of file containing regions to include,\n"
              << "\t\t                 format: start end\n"
              << "\t\t                 start & end positions should be 1-based inclusive positions.\n"
              << "\t\t--keepGT       : comma separated list of genotype fields\n"
              << "\t\t                 to keep in addition to the GT field.\n"
              << "\t\t--splitMulti   : split multi-allelic sites into multiple bi-allelic sites\n"
              << "\t\t                 Genotypes for samples with the alternate not represented\n"
              << "\t\t                 in each bi-allelic site will be set to '0'.\n"
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
    String gtKeepList = "";
    int minAC = -1;
    bool allfields = false;
    bool uncompress = false;
    bool splitMulti = false;
    bool params = false;

    IntervalTree<int> regions;
    std::vector<int> intersection;

    // Read in the parameters.    
    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inputVcf)
        LONG_STRINGPARAMETER("out", &outputVcf)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_PARAMETER("allfields", &allfields)
        LONG_PARAMETER("uncompress", &uncompress)
        LONG_STRINGPARAMETER("sampleSubset", &sampleSubset)
        LONG_INTPARAMETER("minAC", &minAC)
        LONG_STRINGPARAMETER("filterList", &filterList)
        LONG_STRINGPARAMETER("keepGT", &gtKeepList)
        LONG_PARAMETER("splitMulti", &splitMulti)
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
    
    if(!allfields)
    {
        // Do not keep any records with a missing GT or with a non-phased sample.
        inFile.addDiscardRules(VcfFileReader::DISCARD_NON_PHASED | 
                               VcfFileReader::DISCARD_MISSING_GT);
    }

    if(!filterList.IsEmpty())
    {
        // Open the filter list.
        IFILE regionFile = ifopen(filterList, "r");
        String regionLine;
        StringArray regionColumn;
        int start;
        int end;
        int intervalVal = 1;
        if(regionFile == NULL)
        {
            std::cerr << "Failed to open " << filterList 
                      << ", so keeping all positions\n";
            filterList.Clear();
        }
        else
        {
            while( regionFile->isOpen() && !regionFile->ifeof())
            {
                // Read the next interval
                regionLine.Clear();
                regionLine.ReadLine(regionFile);
                if(regionLine.IsEmpty())
                {
                    // Nothing on this line, continue to the next.
                    continue;
                }
                regionColumn.ReplaceColumns(regionLine, ' ');
                if(regionColumn.Length() != 2)
                {
                    std::cerr << "Improperly formatted region line: " 
                              << regionLine << "; skipping to the next line.\n";
                    continue;
                }
                // Convert the columns to integers.
                if(!regionColumn[0].AsInteger(start))
                {
                    // The start position (1st column) is not an integer.
                    std::cerr << "Improperly formatted region line, start position "
                              << "(1st column) is not an integer: "
                              << regionColumn[0]
                              << "; Skipping to the next line.\n";
                    continue;
                }
                if(!regionColumn[1].AsInteger(end))
                {
                    // The start position (1st column) is not an integer.
                    std::cerr << "Improperly formatted region line, end position "
                              << "(2nd column) is not an integer: "
                              << regionColumn[1]
                              << "; Skipping to the next line.\n";
                    continue;
                }
                // Add 1-based inclusive intervals.
                regions.add(start,end, intervalVal);
            }
        }
    }

    int numReadRecords = 0;
    int numWrittenRecords = 0;
    int returnVal = 0;

    std::cerr << "Starting VCF reference panel generation... \n";

    // If not all fields are kept, only keep the specified ones.
    if(!allfields)
    {
        // Set to only store/write the GT field.
        VcfRecordGenotype::addStoreField("GT");

        // Parse the keep list to store/write those fields too.
        if(!gtKeepList.IsEmpty())
        {
            // Keep additional fields.
            StringArray gtFields;
            gtFields.ReplaceColumns(gtKeepList, ',');
            for(int i = 0; i < gtFields.Length(); i++)
            {
                VcfRecordGenotype::addStoreField(gtFields[i].c_str());
            }
        }
    }

    ReusableVector<std::string> origAltArray;
    std::vector<std::vector<std::pair<int, int> > > gtVals;
    std::pair<int, int> gtPair;

    while(inFile.readRecord(record))
    {
        if(!filterList.IsEmpty())
        {
            // Check if the region should be kept.
            intersection.clear();
            regions.get_intersecting_intervals(record.get1BasedPosition(), intersection);
            
            if(intersection.empty())
            {
                // not in the interval, so continue to the next record.
                continue;
            }
        }

        ++numReadRecords;
        
        if(!allfields)
        {
            // Clear the INFO field if not all fields are kept.
            record.getInfo().clear();
        }

        if(splitMulti && (record.getNumAlts() > 1))
        {
            origAltArray.clear();
            gtVals.clear();
            gtVals.resize(record.getNumAlts());
            for(unsigned int i = 1; i <= record.getNumAlts(); i++)
            {
                // Alts start at index 1 of getAlleles
                origAltArray.getNextEmpty() = record.getAlleles(i);
                
            }
            // Loop through samples to store orig GT
            for(int smNum = 0; smNum < record.getNumSamples(); smNum++)
            {
                // For each sample, loop through the GTs.
                for(int gtNum = 0; gtNum < record.getNumGTs(smNum); gtNum++)
                {
                    int allele = record.getGT(smNum, gtNum);
                    if(allele > 0)
                    {
                        // Alternate, store it.
                        gtPair.first = smNum;
                        gtPair.second = gtNum;
                        gtVals[allele-1].push_back(gtPair);
                        // Init all to 0
                        record.setGT(smNum, gtNum, 0);
                    }
                }
            }
            // Need to copy to string since const char* points to the record
            // value which changes when a new ID is set.
            std::string origID = record.getIDStr();
            std::string newID;
            // Loop through each alt (start at 0 of origAltArray).
            for(int i = 0; i < origAltArray.size(); i++)
            {
                record.setAlt(origAltArray.get(i).c_str());
                newID = origID;
                newID += '_' + origAltArray.get(i);
                record.setID(newID.c_str());
                // Loop through and update GTs for this alt.
                for(unsigned int j = 0; j < gtVals[i].size(); j++)
                {
                    record.setGT(gtVals[i][j].first, gtVals[i][j].second, 1);
                }

                // Write the record.
                if(!outFile.writeRecord(record))
                {
                    // Write error.
                    std::cerr << "Failed writing a vcf record.\n";
                    returnVal = -1;
                }
                // Loop through and reset GTs for this alt back to 0.
                for(unsigned int j = 0; j < gtVals[i].size(); j++)
                {
                    record.setGT(gtVals[i][j].first, gtVals[i][j].second, 0);
                }
                ++numWrittenRecords;
            }
            
        }
        else
        {
            // Write the record.
            if(!outFile.writeRecord(record))
            {
                // Write error.
                std::cerr << "Failed writing a vcf record.\n";
                returnVal = -1;
            }
            ++numWrittenRecords;
        }
    }
 
    std::cerr << "NumReadRecords: " << inFile.getNumRecords()
              << "; NumRecordsDiscarded: "
              << inFile.getNumRecords() - numReadRecords
              << "; NumWrittenRecords: " << numWrittenRecords << "\n";

    inFile.close();   
    outFile.close();   

    return(returnVal);
}
