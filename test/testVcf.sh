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


../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileRegionAll.vcf --uncompress --allfield --filterList testFiles/regionList.txt 2> results/vcfFileRegionAll.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileGQHQall.vcf --keepGT GQ,HQ --allfield --uncompress 2> results/vcfFileGQHQall.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileSubset1MinAC2all.vcf --allfield --uncompress --minAC 2 --sampleSubset testFiles/subset1.txt 2> results/vcfFileSubset1MinAC2all.vcf.log

../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileBi.vcf --uncompress --splitMulti 2> results/vcfFileBi.vcf.log
../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileSubset1allSplit.vcf --allfield --uncompress --sampleSubset testFiles/subset1.txt --splitMulti 2> results/vcfFileSubset1allSplit.vcf.log

../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileSubset1MinAC1allSplit.vcf --minAC 1 --allfield --uncompress --sampleSubset testFiles/subset1.txt --splitMulti 2> results/vcfFileSubset1MinAC1allSplit.vcf.log

../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileSubset1MinAC2allSplit.vcf --minAC 2 --allfield --uncompress --sampleSubset testFiles/subset1.txt --splitMulti 2> results/vcfFileSubset1MinAC2allSplit.vcf.log

../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileSubset1MinAC2allSplitUpdateID.vcf --minAC 2 --allfield --uncompress --sampleSubset testFiles/subset1.txt --splitMulti --idUpdate 2> results/vcfFileSubset1MinAC2allSplitUpdateID.vcf.log

../bin/vcfRefGen --in testFiles/vcfFile.vcf --out results/vcfFileBiUpdateID.vcf --uncompress --splitMulti --idUpdate 2> results/vcfFileBiUpdateID.vcf.log

../bin/vcfRefGen --in testFiles/vcfFileSV.vcf --out results/vcfFileSVBiUpdateID.vcf --uncompress --splitMulti --idUpdate 2> results/vcfFileSVBiUpdateID.vcf.log

diff expected results