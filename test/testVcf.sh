../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFile.vcf --uncompress 2> results/vcfFile.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFile.vcf.gz 2> results/vcfFile.vcf.gz.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileSubset1.vcf --uncompress --sampleSubset testFiles/subset1.txt 2> results/vcfFileSubset1.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileMinAC2.vcf --uncompress --minAC 2 2> results/vcfFileMinAC2.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileSubset1MinAC2.vcf --uncompress --minAC 2 --sampleSubset testFiles/subset1.txt 2> results/vcfFileSubset1MinAC2.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileMinAC1.vcf --uncompress --minAC 1 2> results/vcfFileMinAC1.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileSubset1MinAC1.vcf --uncompress --minAC 1 --sampleSubset testFiles/subset1.txt 2> results/vcfFileSubset1MinAC1.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileRegion.vcf --uncompress --filterList testFiles/regionList.txt 2> results/vcfFileRegion.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileDP.vcf --keepGT DP --uncompress 2> results/vcfFileDP.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileGQHQ.vcf --keepGT GQ,HQ --uncompress 2> results/vcfFileGQHQ.vcf.log



diff expected results