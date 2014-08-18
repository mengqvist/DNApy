#!/usr/bin/python
# Copyright (C) 2014 Brian J. Stucky
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Script by Brian J. Stucky was slightly modified by Martin KM Engqvist

import SequenceTrace


class TraceFileError(Exception):
    pass

class UnknownFileTypeError(TraceFileError):
    def __str__(self):
        return 'The file format was not recognized.  Please convert the file to a supported format and try again.'


class SCFError(TraceFileError):
    pass

class SCFVersionError(SCFError):
    def __init__(self, version, revision):
        self.version = version
        self.revision = revision

    def __str__(self):
        return 'This file uses version ' + self.version + '.' + self.revision + ' of the SCF format.  This software only supports version 3.00 of the format.'

class SCFDataError(SCFError):
    def __init__(self, expectedlen, actuallen):
        self.expectedlen = expectedlen
        self.actuallen = actuallen

    def __str__(self):
        return 'Error reading SCF file data.  Expected ' + str(self.expectedlen) + ' bytes but only got ' + str(self.actuallen) + ' bytes.  The file appears to be damaged.'


class SCFSequenceTrace(SequenceTrace):
    def loadFile(self, filename):
        self.fname = filename

        try:
            self.tf = open(filename, 'rb')
        except IOError:
            raise

        magicnum = self.tf.read(4)
        #print magicnum
        if magicnum != '.scf':
            raise SCFError('The SCF file header is invalid.  The file appears to be damaged.')

        try:
            numsamps = unpack('>I', self.tf.read(4))[0]
            sampstart = unpack('>I', self.tf.read(4))[0]
            #print numsamps, sampstart

            numbases = unpack('>I', self.tf.read(4))[0]
            # skip 8 bytes
            self.tf.read(8)
            basesstart = unpack('>I', self.tf.read(4))[0]
            #print numbases, basesstart

            commentslen = unpack('>I', self.tf.read(4))[0]
            commentsstart = unpack('>I', self.tf.read(4))[0]
            #print commentslen, commentsstart

            version = self.tf.read(4)
            #print version

            samplesize = unpack('>I', self.tf.read(4))[0]
            codeset = unpack('>I', self.tf.read(4))[0]
            #print samplesize, codeset
        except struct.error:
            raise SCFError('The SCF file header is invalid.  The file appears to be damaged.')

        if version != '3.00':
            raise SCFVersionError(version[0], version[2:])

        if samplesize not in (1, 2):
            raise SCFError('Invalid sample size value in SCF header.  The size specified was ' + str(samplesize) + ', but must be either 1 or 2.')

        if codeset != 0:
            raise SCFError('Invalid code set specified in SCF header.  This file uses code set ' + str(codeset) + ', but this software only supports code set 0.')

        self.readBasesData(numbases, basesstart)
        self.readTraceData(numsamps, sampstart, samplesize)
        self.readComments(commentslen, commentsstart)

    def readBasesData(self, numbases, basesstart):
        self.basepos = list()
        probs = {'A': list(), 'C': list(), 'G': list(), 'T': list()}
        self.bcconf = list()

        self.tf.seek(basesstart, 0)

        try:
            # get the base locations
            for cnt in range(0, numbases):
                index = unpack('>I', self.tf.read(4))[0]
                self.basepos.append(index)
            #print self.basepos

            # get the base call probabilities for all bases
            for base in ('A', 'C', 'G', 'T'):
                for cnt in range(0, numbases):
                    prob = unpack('B', self.tf.read(1))[0]
                    probs[base].append(prob)
        except struct.error:
            raise SCFError('Error while reading base call locations and probabilities from the SCF file.  The file appears to be damaged.')

        # get the base calls
        self.basecalls = self.tf.read(numbases).upper()
        #print self.basecalls
        if numbases != len(self.basecalls):
            raise SCFDataError(numbases, len(self.basecalls))

        # build the confidence scores list
        for cnt in range(0, numbases):
            base = self.basecalls[cnt]
            self.bcconf.append(probs[base][cnt])

        #print self.bcconf

    def readTraceData(self, numsamps, sampstart, sampsize):
        if sampsize == 1:
            formatstr = 'B'
        else:
            formatstr = '>H'

        self.tf.seek(sampstart, 0)

        maxval = 0

        for base in ('A', 'C', 'G', 'T'):
            samps = list()

            try:
                # read the raw sample data
                for cnt in range(0, numsamps):
                    val = unpack(formatstr, self.tf.read(sampsize))[0]
                    samps.append(val)
            except struct.error:
                raise SCFDataError((numsamps * sampsize), (len(samps) * sampsize))

            # sample values are double-delta encoded (i.e., two successive rounds of differences)
            if sampsize == 1:
                self.decode8BitDoubleDelta(samps)
            else:
                self.decode16BitDoubleDelta(samps)

            self.tracesamps[base] = samps
            tmpmax = max(self.tracesamps[base])
            if tmpmax > maxval:
                maxval = tmpmax

        self.max_traceval = maxval
        #print self.tracesamps['A']

    def decode8BitDoubleDelta(self, data):
        for clev in range(0, 2):
            prev = 0
            for cnt in range(0, len(data)):
                actual = data[cnt] + prev

                # simulate 1-byte integer overflow, if needed
                if actual > 255:
                    actual -= 256

                prev = actual
                data[cnt] = actual
    
    def decode16BitDoubleDelta(self, data):
        for clev in range(0, 2):
            prev = 0
            for cnt in range(0, len(data)):
                actual = data[cnt] + prev

                # simulate 2-byte integer overflow, if needed
                if actual > 65535:
                    actual -= 65536

                prev = actual
                data[cnt] = actual

    def readComments(self, commentslen, commentsstart):
        self.tf.seek(commentsstart, 0)

        total = 0
        while total < (commentslen - 1):
            line = self.tf.readline()
            total += len(line)

            if line == '':
                raise SCFError('Unable to read the comments section of the SCF file.  The file appears to be damaged.')

            # get rid of the trailing '\n'
            line = line[:-1]

            key, sep, value = line.partition('=')
            #print key + ': ' + value
            self.comments[key] = value

        # make sure the next character is the null-terminator for the comments list
        if self.tf.read(1) != '\0':
            raise SCFError('Missing null character at end of comments section.  The file appears to be damaged.')




