﻿using Bio;
using Bio.IO.FastA;
using NUnit.Framework;
using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;
using WorkflowLayer;

namespace Test
{
    [TestFixture]
    public class WrapperIntegrationTests
    {
        #region Installs

        [Test, Order(0)]
        public void TestInstall()
        {
            ManageToolsFlow.Install(TestContext.CurrentContext.TestDirectory);

            // bedops
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "bedops")));

            // bedtools
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "bedtools2", "bin", "bedtools")));

            // cufflinks
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "cufflinks-2.2.1")));

            // gatk
            Assert.IsTrue(Directory.GetFiles(Path.Combine(TestContext.CurrentContext.TestDirectory, "gatk"), "gatk*local.jar").Length > 0);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "ChromosomeMappings")));

            // hisat2
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "hisat2-2.1.0", "hisat2")));

            // mfold
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "mfold-3.6", "scripts", "mfold")));

            // rsem
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSEM-1.3.0", "rsem-prepare-reference")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSEM-1.3.0", "rsem-calculate-expression")));

            // rseqc
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSeQC-2.6.4")));

            // samtools
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "samtools-1.6", "samtools")));

            // scalpel
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "scalpel-0.5.3")));

            // skewer
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "skewer-0.2.2")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "BBMap")));

            // slncky
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "slncky")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "slncky", "annotations")));
            Assert.IsTrue(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "lastz*").Length > 0);

            // snpeff
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "snpEff", "snpEff.jar")));

            // sratoolkit
            Assert.IsTrue(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "sratoolkit*").Length > 0);
            Assert.IsTrue(Directory.GetFiles(Directory.GetDirectories(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "sratoolkit*")[0], "bin")[0], "fastq-dump").Length > 0);

            // star
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR")));

            // star-fusion
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR-Fusion_v1.1.0")));
        }

        private string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa");

        [Test, Order(1)]
        public void DownloadReferences()
        {
            EnsemblDownloadsWrapper downloadsWrapper = new EnsemblDownloadsWrapper();
            downloadsWrapper.DownloadReferences(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37");

            // - a basic set of chromosomes, fairly small ones
            string a = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            // - chromosomes and contigs that test ordering: 9 comes before 22 in karyotipic order, but not lexographic
            string b = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.fa");

            if (!File.Exists(a) || !File.Exists(b))
            {
                List<ISequence> chromosomes = new FastAParser().Parse(new FileStream(genomeFastaPath, FileMode.Open)).ToList();
                FastAFormatter formatter = new FastAFormatter();

                if (!File.Exists(a))
                    Genome.WriteFasta(chromosomes.Where(x => x.ID.StartsWith("20") || x.ID.StartsWith("21") || x.ID.StartsWith("22")), a);

                if (!File.Exists(b))
                    Genome.WriteFasta(chromosomes.Where(x => x.ID.StartsWith("9") || x.ID.StartsWith("22") || x.ID.StartsWith("GL000210") || x.ID.StartsWith("HG1287_PATCH")), b);
            }

            // Additional setup for small integration tests
            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "setup.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                "if [ -f Homo_sapiens.GRCh37.75.gtf.gz ]; then gunzip Homo_sapiens.GRCh37.75.gtf.gz; fi",
                @"if [ ! -f 202122.gtf ]; then grep '^20\|^21\|^22' Homo_sapiens.GRCh37.75.gtf > 202122.gtf; fi",
                @"if [ ! -f 922HG1287_PATCH.gtf ]; then grep '^9\|^22\|^HG1287_PATCH\|^GL000210.1' Homo_sapiens.GRCh37.75.gtf > 922HG1287_PATCH.gtf; fi",
            }).WaitForExit();

            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.fa")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.gtf")));
        }

        #endregion Installs

        #region SRA download test

        [Test, Order(1)]
        public void TestDownloadSRA()
        {
            SRAToolkitWrapper sratoolkit = new SRAToolkitWrapper();
            sratoolkit.Fetch(TestContext.CurrentContext.TestDirectory, "SRR6304532", TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(sratoolkit.FastqPaths.All(f => File.Exists(f)));
            Assert.IsTrue(File.Exists(sratoolkit.LogPath));
        }

        #endregion SRA download test

        #region BED conversion tests

        [Test, Order(1)]
        public void TestConvertGff()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(1)]
        public void TestConvertGtf()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(1)]
        public void TestConvertGtf12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed12")).Length > 0);
        }

        [Test, Order(1)]
        public void TestConvertGffToBed12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gff.gff3") + ".bed12")).Length > 0);
        }

        #endregion BED conversion tests

        

        #region RSEM Tests

        [Test, Order(2)]
        public void RSEMStarCalculate()
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperRsemStar.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                RSEMAlignerOption.STAR,
                Strandedness.None,
                new[] { newMapper },
                true);
            quantification.QuantifyTranscripts();

            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GeneResultsSuffix));
            Assert.IsTrue(Directory.Exists(quantification.RsemOutputPrefix + RSEMWrapper.StatDirectorySuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TimeSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TranscriptBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamIndexSuffix));
        }

        [Test, Order(2)]
        public void RSEMBowtieCalculate()
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperRsemBowtie.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                RSEMAlignerOption.Bowtie2,
                Strandedness.None,
                new[] { newMapper },
                true);
            quantification.QuantifyTranscripts();

            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GeneResultsSuffix));
            Assert.IsTrue(Directory.Exists(quantification.RsemOutputPrefix + RSEMWrapper.StatDirectorySuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TimeSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TranscriptBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamIndexSuffix));
        }

        [Test, Order(4)]
        public void RSEMStarCalculateFromPaired()
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperRsemBowtie.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                RSEMAlignerOption.STAR,
                Strandedness.None,
                new[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SRR6319804_1-trimmed-pair1.segment.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SRR6319804_1-trimmed-pair2.segment.fastq"),
                },
                true);
            quantification.QuantifyTranscripts();

            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GeneResultsSuffix));
            Assert.IsTrue(Directory.Exists(quantification.RsemOutputPrefix + RSEMWrapper.StatDirectorySuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TimeSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TranscriptBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamIndexSuffix));
        }

        // I'm having trouble getting RSEM to work with comma-separated inputs... I think it's because of STAR, which I have had trouble with in this respect in the past.

        //[Test, Order(2)]
        //public void RSEMStarCalculate2Fastq()
        //{
        //    string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperRsemStar2.fastq");
        //    if (!File.Exists(newMapper))
        //    {
        //        File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper.fastq"), newMapper);
        //    }
        //    string newMapper2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperRsemStar2Again.fastq");
        //    if (!File.Exists(newMapper2))
        //    {
        //        File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper.fastq"), newMapper2);
        //    }

        //    TranscriptQuantificationFlow.QuantifyTranscripts(
        //        TestContext.CurrentContext.TestDirectory,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
        //        Environment.ProcessorCount,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
        //        RSEMAlignerOption.STAR,
        //        Strandedness.None,
        //        new[] { newMapper + "," + newMapper2 },
        //        true,
        //        out string referencePrefix,
        //        out string outputPrefix);

        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.IsoformResultsSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GeneResultsSuffix));
        //    Assert.IsTrue(Directory.Exists(outputPrefix + RSEMWrapper.StatDirectorySuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TimeSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TranscriptBamSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TranscriptBamIndexSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GenomeBamSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GenomeBamIndexSuffix));
        //}

        //[Test, Order(4)]
        //public void RSEMStarCalculateTwoPairFastq()
        //{
        //    TranscriptQuantificationFlow.QuantifyTranscripts(
        //        TestContext.CurrentContext.TestDirectory,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
        //        Environment.ProcessorCount,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
        //        RSEMAlignerOption.STAR,
        //        Strandedness.None,
        //        new[]
        //        {
        //            String.Join(",", new string[] {
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_1.fastq"),
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_1.fastq") }),
        //            String.Join(",", new string[] {
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_2.fastq"),
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_2.fastq") })
        //        },
        //        true,
        //        out string referencePrefix,
        //        out string outputPrefix);

        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.IsoformResultsSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GeneResultsSuffix));
        //    Assert.IsTrue(Directory.Exists(outputPrefix + RSEMWrapper.StatDirectorySuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TimeSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TranscriptBamSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TranscriptBamIndexSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GenomeBamSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GenomeBamIndexSuffix));
        //}

        #endregion RSEM Tests

        #region Infer Experiment tests

        [Test, Order(4)]
        public void StrandSpecificityTest()
        {
            BAMProperties bam = new BAMProperties(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperAgain-trimmedAligned.sortedByCoord.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")),
                0.8);
            Assert.AreEqual(Strandedness.None, bam.Strandedness);
            Assert.AreEqual(RnaSeqProtocol.SingleEnd, bam.Protocol);
        }

        [Test, Order(2)]
        public void InnerDistanceTest()
        {
            Assert.AreEqual(132, RSeQCWrapper.InferInnerDistance(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "paired_end.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                out string[] outputFiles));
        }

        #endregion Infer Experiment tests

        #region Skewer tests

        [Test, Order(2)]
        public void SkewerSingle()
        {
            SkewerWrapper.Trim(TestContext.CurrentContext.TestDirectory,
                1,
                19,
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "read1.fastq") },
                out string[] readTrimmedPaths,
                out string log);
            Assert.True(File.Exists(readTrimmedPaths[0]));
            Assert.True(File.Exists(log));
            File.Delete(readTrimmedPaths[0]);
            File.Delete(log);
        }

        [Test, Order(2)]
        public void SkewerPaired()
        {
            SkewerWrapper.Trim(TestContext.CurrentContext.TestDirectory,
                19,
                1,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","read2.fastq")
                },
                out string[] readTrimmedPaths,
                out string log);
            Assert.True(File.Exists(readTrimmedPaths[0]));
            Assert.True(File.Exists(readTrimmedPaths[1]));
            Assert.True(File.Exists(log));
            File.Delete(readTrimmedPaths[0]);
            File.Delete(readTrimmedPaths[1]);
            File.Delete(log);
        }

        [Test, Order(2)]
        public void SkewerPairedGz()
        {
            SkewerWrapper.Trim(TestContext.CurrentContext.TestDirectory,
                19,
                1,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "read1.fastq.gz"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","read2.fastq.gz")
                },
                out string[] readTrimmedPaths,
                out string log);
            Assert.True(Path.GetFileName(readTrimmedPaths[0]) == "read1-trimmed-pair1.fastq");
            Assert.True(Path.GetFileName(readTrimmedPaths[1]) == "read1-trimmed-pair2.fastq");
            Assert.True(Path.GetFileName(log) == "read1-trimmed.log");
            File.Delete(readTrimmedPaths[0]);
            File.Delete(readTrimmedPaths[1]);
            File.Delete(log);
        }

        #endregion Skewer tests

        #region GATK tests

        [Test, Order(2)]
        public void DownloadKnownSites()
        {
            GATKWrapper.DownloadEnsemblKnownVariantSites(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                true,
                "grch37",
                out string ensemblKnownSitesPath);
            Assert.IsTrue(File.Exists(ensemblKnownSitesPath));

            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "setupKnownSitesTest.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf")) + " ]; " +
                    @"then grep '^#\|^chr20\|^chr21\|^chr22\|^20\|^21\|^22' " + WrapperUtility.ConvertWindowsPath(ensemblKnownSitesPath) +
                    " > 202122.vcf; fi",
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf")) + " ]; " +
                    @"then grep '^#\|^chr9\|^chr22\|^chrHG1287_PATCH\|chr21_gl000210_random\|^9\|^22\|^HG1287_PATCH\|^GL000210.1' " + WrapperUtility.ConvertWindowsPath(ensemblKnownSitesPath) +
                    " > 922HG1287_PATCH.vcf; fi",
            }).WaitForExit();
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf")));
        }

        [Test, Order(4)]
        public void GatkWorflow()
        {
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "gatkWorkflowTest.bash"), new List<string>
            (
                GATKWrapper.PrepareBamAndFasta(TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperAgain-trimmedAligned.sortedByCoord.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                "grch37",
                out string new_bam).Concat(

                GATKWrapper.SplitNCigarReads(TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                new_bam,
                out string splitTrimBam).Concat(

                // No longer needed with HaplotypeCaller
                //GATKWrapper.RealignIndels(TestContext.CurrentContext.TestDirectory,
                //    8,
                //    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                //    new_bam,
                //    out string realigned_bam,
                //    ""); // not including known sites speeds this up substantially, and I'm not planning to use these indels

                // Takes kind of a long time, and it's not recommended for RNA-Seq yet
                //GATKWrapper.base_recalibration(TestContext.CurrentContext.TestDirectory,
                //    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "202122.fa"),
                //    realigned_bam,
                //    out string recal_table_filepath,
                //    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "202122.vcf"));

                GATKWrapper.VariantCalling(TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                splitTrimBam,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf"),
                out string newVcf)))
            )).WaitForExit();

            Assert.IsTrue(File.Exists(newVcf));
            Assert.IsTrue(new FileInfo(newVcf).Length > 0);
        }

        [Test, Order(5)]
        public void convertVcf()
        {
            GATKWrapper.ConvertVCFChromosomesUCSC2Ensembl(TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1Ucsc.vcf"),
                "grch37",
                out string ensemblVcf);
            Assert.IsTrue(File.Exists(ensemblVcf));
            Assert.IsTrue(new FileInfo(ensemblVcf).Length > 0);
        }

        #endregion GATK tests

        #region Cufflinks and Stringtie tests

        [Test, Order(4)]
        public void CufflinksRun()
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.sortedByCoord.out.bam");
            string script_name = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "cufflinksRun.bash");
            WrapperUtility.GenerateAndRunScript(script_name, CufflinksWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")),
                false,
                false,
                out string outputDirectory
                )).WaitForExit();
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, CufflinksWrapper.TranscriptsFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, CufflinksWrapper.SkippedTranscriptsFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, CufflinksWrapper.IsoformAbundanceFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, CufflinksWrapper.GeneAbundanceFilename)));
        }

        [Test, Order(4)]
        public void StringtieRun()
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.sortedByCoord.out.bam");
            string script_name = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "stringtieRun.bash");
            WrapperUtility.GenerateAndRunScript(script_name, StringTieWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")),
                Strandedness.None,
                false,
                out string outputGtf
                )).WaitForExit();
            Assert.IsTrue(File.Exists(outputGtf));
            Assert.IsTrue(new FileInfo(outputGtf).Length > 0);
        }

        #endregion Cufflinks tests

        #region Slncky tests

        [Test, Order(5)]
        public void SlnckyRun()
        {
            string cufflinksTranscripts = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.sortedByCoord.out.cufflinksOutput", CufflinksWrapper.TranscriptsFilename);
            string slnckyOutPrefix = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "cuffmergeModel848835768.slnckyOut", "annotated"); // strange folder, so that it covers the same test as the lncRNAdiscovery run
            string scriptName = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "SlnckyRun.bash");
            WrapperUtility.GenerateAndRunScript(scriptName, 
                SlnckyWrapper.Annotate(
                    TestContext.CurrentContext.TestDirectory,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                    Environment.ProcessorCount,
                    cufflinksTranscripts,
                    "GRCh37",
                    slnckyOutPrefix
                )).WaitForExit();
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.CanonicalToLncsSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.ClusterInfoSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.FilteredInfoSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.LncsBedSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.LncsInfoSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrfsSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrthologsSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrthologsTopSuffix));
        }

        #endregion Slncky tests

        #region Scalpel tests

        [Test, Order(4)]
        public void ScalpelCall()
        {
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "scalpel.bash"),  
                ScalpelWrapper.CallIndels(TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.bed12"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.sortedByCoord.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "scalpel_test_out"),
                out string newVcf)).WaitForExit();
            Assert.IsTrue(File.Exists(newVcf));
        }

        #endregion Scalpel tests

        #region SnpEff tests

        [Test, Order(1)]
        public void DownloadSnpEffDatabase()
        {
            SnpEffWrapper.DownloadSnpEffDatabase(TestContext.CurrentContext.TestDirectory,
                "grch37",
                out string databaseListPath);
            Assert.IsTrue(File.Exists(databaseListPath));
            string[] databases = Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "snpEff", "data"));
            Assert.IsTrue(databases.Any(x => Path.GetFileName(x).StartsWith("grch37", true, null)));
        }

        [Test, Order(4)]
        public void BasicSnpEffAnnotation()
        {
            SnpEffWrapper.PrimaryVariantAnnotation(TestContext.CurrentContext.TestDirectory,
                "grch37",
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.sortedByCoord.outProcessed.out.fixedQuals.split.vcf"),
                out string html,
                out string annVcf,
                out string annGenes
                );
            Assert.IsTrue(File.Exists(html) && new FileInfo(html).Length > 0);
            Assert.IsTrue(File.Exists(annVcf) && new FileInfo(annVcf).Length > 0);
            Assert.IsTrue(File.Exists(annVcf) && new FileInfo(annGenes).Length > 0);
        }

        #endregion SnpEff tests

        #region Runner Tests

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single-end
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromFastqs()
        {
            SampleSpecificProteinDBFlow ssdbf = new SampleSpecificProteinDBFlow();
            ssdbf.Parameters = new SampleSpecificProteinDBParameters(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                Environment.ProcessorCount,
                new List<string[]>
                {
                    new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper.fastq") },
                    new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperAgain.fastq") },
                },
                false,
                false,
                false,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf"));
            ssdbf.GenerateSAVProteinsFromFastqs();
            foreach (string database in ssdbf.VariantAppliedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains(FunctionalClass.MISSENSE.ToString())));
            }
        }

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single-end
        /// </summary>
        [Test, Order(4)]
        public void LncRnaDiscoveryRunFromFastqs()
        {
            LncRNADiscoveryFlow lncRNAdiscovery = new LncRNADiscoveryFlow();
            lncRNAdiscovery.Parameters = new LncRNADiscoveryParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37",
                Environment.ProcessorCount,
                new List<string[]>
                {
                    new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper.fastq") },
                },
                false,
                false,
                false,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                true);
            lncRNAdiscovery.LncRNADiscoveryFromFastqs();

            Assert.IsTrue(File.Exists(lncRNAdiscovery.MergedGtfPath));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.CanonicalToLncsSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.ClusterInfoSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.FilteredInfoSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.LncsBedSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.LncsInfoSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.OrfsSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.OrthologsSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.OrthologsTopSuffix));
        }

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single end
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromTwoPairsFastqs()
        {
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_1.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_1.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_1.fastq"));
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_2.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_2.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_2.fastq"));

            SampleSpecificProteinDBFlow ssdbf = new SampleSpecificProteinDBFlow();
            ssdbf.Parameters = new SampleSpecificProteinDBParameters(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                Environment.ProcessorCount,
                new List<string[]>
                {
                    new string[] {
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_1.fastq"),
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_2.fastq") },
                    new string[] {
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_1.fastq"),
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_2.fastq") }
                },
                false,
                false,
                false,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf"));

            ssdbf.GenerateSAVProteinsFromFastqs();
            foreach (string database in ssdbf.VariantAppliedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                //Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("ANN="))); // no longer see any variations for this test set with variant filtering criteria
            }
        }

        /// <summary>
        /// Handling tough non-karyotypic ordering of chromosomes and an SRA input
        ///
        /// This also tests well-encoded quality scores, so if it starts to fail, check out whether the exit code of the FixMisencodedQualityBaseReads is expected (2 for failure).
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromSRA()
        {
            var fastqs = SRAToolkitWrapper.GetFastqsFromSras(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), "SRR6319804");
            SampleSpecificProteinDBFlow ssdbf = new SampleSpecificProteinDBFlow();
            ssdbf.Parameters = new SampleSpecificProteinDBParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37",
                Environment.ProcessorCount,
                fastqs,
                false,
                false,
                false,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.gtf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf"), // there is no equivalent of the patch; just checking that that works
                null,
                0.05,
                7,
                true,
                5000);
            ssdbf.GenerateSAVProteinsFromFastqs();
                
            foreach (string database in ssdbf.VariantAppliedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("variant")));// no variants anymore with the filtering criteria
            }
        }

        #endregion Runner Tests
    }
}