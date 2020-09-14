"""
ePSproc ePSdata downloader

Basic IO function for downloading data from ePSdata, via Zenodo.

20/07/20    v1, basics in place.

TODO:

* More integration with ePSman. Use of json info files?
* Multipart zip file handling.
* Download multiple jobs/files. Could be in class, or as object per Zenodo record. Will be useful for, e.g., cases comparing multiple orbitals.
* Tidy up.

"""

#**** Imports
# Core
import requests
import wget
import zipfile
import os
import glob
# import shutil
from pathlib import Path

from collections import Counter
import pprint

# Web stuff
from urllib.parse import urlparse
from IPython.core.display import HTML  # For HTML rendering support
from html.parser import HTMLParser

#**** Local util functions
# TODO: move these elsewhere, for general use.

# Basic bytes to KB/Mb... conversion, from https://stackoverflow.com/questions/2104080/how-to-check-file-size-in-python
def convert_bytes(num):
    """
    This function will convert bytes to MiB, GiB... etc

    From https://stackoverflow.com/a/39988702

    Thanks to @rajiv-sharma: https://stackoverflow.com/users/2679465/rajiv-sharma

    Converted to binary prefix units (https://superuser.com/questions/1076888/looking-for-clarification-on-binary-prefix-logic-history-vs-si-prefix/1077275#1077275)

    """

    for x in ['bytes', 'KiB', 'MiB', 'GiB', 'TiB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
#             return [num, x]
        num /= 1024.0


# Basic URL from HTML parser, from https://stackoverflow.com/a/6883228
class hrefParser(HTMLParser):
    """
    Basic URL from HTML parser class, from https://stackoverflow.com/a/6883228

    Thanks to @senderle: https://stackoverflow.com/users/577088/senderle

    Example
    -------

    >>> p = hrefParser()
    >>> p.feed(myString)
    >>> p.output_list

    """
    def __init__(self, output_list=None):
        HTMLParser.__init__(self)
        if output_list is None:
            self.output_list = []
        else:
            self.output_list = output_list
    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.output_list.append(dict(attrs).get('href'))


# General file unzip handling, including multi-part zips
def unzipFile(fileIn, unzipPath = None, confirmFlag = True, verbose = True):
    """
    Basic zip file handling

    Parameters
    ----------
    fileIn : str or Path object
        Zip file to unzip.

    unzipPath : str or Path object, optional, default = None
        Path to unzip files to.
        If not passed, will be set as fileIn directory.

    TODO:
    Add file size checking and confirmation in case of large files.

    """
    if unzipPath is None:
        unzipPath = Path(fileIn).parent

    with zipfile.ZipFile(fileIn,"r") as zipObj:

        # Get details
        zipFiles = zipObj.namelist()
        zipInfo = zipObj.infolist()
        zipSize = sum([item.file_size for item in zipInfo])

        if verbose:
            print(f'\n*** Unzipping archive: {fileIn}')
            print(f'Unzipped archive size will be {convert_bytes(zipSize)}.')

        # Manually confirm unzip if option set.
        unzipFlag = True
        if confirmFlag:
            if input('Unzip? (y/n): ') == 'y':
                unzipFlag = True
            else:
                unzipFlag = False

        if unzipFlag:
            zipObj.extractall(unzipPath)

            if verbose:
                print(f"Unzipped file {fileIn} to directory {unzipPath}")
        else:
            if verbose:
                print('Skipped unzipping.')

        zipDetails={'path':unzipPath,
                    'zipfile':fileIn,
                    'files':zipFiles,
                    'info':zipInfo,
                    'unzipped':unzipFlag}

    return zipDetails


# Basic handling for multi-part zip files.
def joinZip(fileIn, verbose = True):
    """
    Join multipart zip file from parts [.zip,.z01,.z02....]

    Parameters
    ----------
    fileIn : str
        Zip file to join, pass first in sequence (.zip) only.

    Returns
    -------
    fileOut : str
        Name for joined archive.
        If archive not joined, returns None.

    Note
    ----
    If passed an archive which is not in multiple parts, zip will just copy it.


    TODO:

    - Better os.system usage, return errors etc.
    - Options apart from Linux zip command? Way to do this pure Python?



    """
    # Set output filename (for joined zip)
    fileIn = Path(fileIn)  # For Path object
    fileOut = fileIn.with_name(fileIn.stem + '_joined.zip')

    if verbose:
        print(f'Joining archive parts {fileIn} to {fileOut}.')

    # Check system - only tested on linux using zip command so far.
    if os.uname()[0] == 'Linux':
        zTest = os.system(f"zip -s 0 {fileIn} --out {fileOut}")

        if not zTest:
            return fileOut
        else:
            print(f"*** Multipart zip join for {fileIn} failed.")
            return None

    else:
        print(f"*** Multipart zip join not implemented for {os.uname()[0]}")

        return None

    # Basics from epsman code.
        # Just need to run this at command line on remote (Linux)
        # TODO: add some checking logic here, at the moment can fail if fileOut already exists.
        # Note -j for junking the path, see http://manpages.ubuntu.com/manpages/precise/man1/zip.1.html
        # os.system(f"zip -j -s {chunk}m {fileOut} {arch}")

        # TO REBUILD:
        # zip -s 0 testMultipart.zip --out testRecon.zip


# Basic dir prune function to remove extraneous dirs from unzip.
def pruneDirTree(path, verbose = True):
    '''
    Prune empty dirs from dir tree `path` using os.removedirs().

    This will proceed from topmost dir recursively, and stop at first non-empty dir.

    https://docs.python.org/3/library/os.html?highlight=os#os.removedirs

    TODO:
    Might be that another method is preferable for more control/feedback? E.g. pure pathlib recursive method: https://stackoverflow.com/a/49782093

    '''

    try:
        os.removedirs(path)

    # Below is actually redundant for the most part, since os.removedirs() shouldn't return errno 39 (it's ignored).
    # Not sure how to confirm dir pruning in this case.
    # Might be that another method is preferable for more control/feedback.
    except OSError as e:

        # This basically indicates the end of the pruning, so return remaining path.
        if e.errno == 39:
            if verbose:
                print(f'Pruned dir tree back to {e.filename}')

            return e.filename

        else:
            raise



class ePSdata():
    """
    Class for interfacing with ePSdata via Zenodo API.

    Init by passing either DOI, URL or Zenodo ID.

    TODO: add search option here too!  See epsdata for a prototype - may require API key.

    Parameters
    ----------
    doi : str, optional, default=None
        DOI for Zenodo record, e.g. '10.5281/zenodo.3629721'

    URL : str, optional, default=None
        URL for Zenodo record, can be DOI or Zenodo address,
        e.g. 'http://dx.doi.org/10.5281/zenodo.3629721' or 'https://zenodo.org/record/3629721'

    ID : str or int, optional, default=None
        Zenodo ID number for record, e.g. 3629721

    downloadDir : str or Path object, optional, default=None
        Path for downloaded file(s).
        Defaults to current working dir, from os.getcwd().

    """

    def __init__(self, doi=None, ID=None, URL=None, downloadDir=None):

        # Set record
        self.setRecord(doi=doi, ID=ID, URL=URL)
        self.setLocal(downloadBase=downloadDir)

        # Get record info
        self.getInfo()

        # Make download dir if remote record OK
        if self.r.ok:
            self.makeLocal()


    def setRecord(self, doi=None, ID=None, URL=None):
        # Set base locations - may want to do this elsewhere, and/or with looping dicts?
        doiBase = '10.5281/zenodo.'  # Common Zenodo DOI base
        doiURLBase = 'http://dx.doi.org/'
        zenURLBase = 'https://zenodo.org/record/'  # zenURLBase + ID corresponds to record webpage.
        zenAPIBase = 'https://zenodo.org/api/records/'  # Need /api for API calls!

        self.base = {}  # TODO - set above to dict item.

        # Set record IDs, starting from DOI
        # TODO - tidy & check logic/repetition here. Bit fugly.
        recordID = {}

        if doi is not None:
            recordID['doi'] = doi  # E.g. '10.5281/zenodo.3629721'
            recordID['url'] = {'doi':doiURLBase + recordID['doi']}
            recordID['zenID'] = int(recordID['doi'].rsplit('.',1)[-1])

        elif URL is not None:
            # Parse URL and check type (doi or Zenodo)
            urlSplit = urlparse(URL)

            if urlSplit.netloc == 'dx.doi.org':
                recordID['url'] = {'doi':URL}
                recordID['doi'] = urlSplit.path.strip('/')  # Set doi from path, remove leading '/'
                recordID['zenID'] = int(urlSplit.path.rsplit('.',1)[-1])


            if urlSplit.netloc == 'zenodo.org':
                recordID['url'] = {'zen':URL}
                recordID['zenID'] = urlSplit.path.rsplit('/')[-1]  # Set ID from URL
                recordID['doi'] = doiBase + recordID['zenID']
                # recordID['url']['doi'] = urljoin(doiURLbase + doiBase, ID)  # With urljoin
                recordID['url']['doi'] = doiURLBase + doiBase + recordID['zenID']  # From strs.


        elif ID is not None:
            recordID['zenID'] = ID
            recordID['doi'] = doiBase + ID
            recordID['url'] = {'zen':zenURLBase + ID,
                               'doi':doiURLBase + doiBase + ID}

        # Set API URL, this will be used to get the record info + files.
        recordID['url']['get'] = zenAPIBase + str(recordID['zenID'])

        self.recordID = recordID


    def setLocal(self, downloadBase=None):
        if downloadBase is None:
            self.recordID['downloadBase'] = Path(os.getcwd()).expanduser()
        else:
            self.recordID['downloadBase'] = Path(downloadBase).expanduser()

        self.recordID['downloadDir'] = self.recordID['downloadBase']/str(self.recordID['zenID'])

        print(f"*** Download dir set to: {self.recordID['downloadDir']}")

    def makeLocal(self):
        try:
            os.mkdir(self.recordID['downloadDir'])
            print(f"\n*** Created {self.recordID['downloadDir']}")
        except FileExistsError:
            print(f"*** Directory {self.recordID['downloadDir']} already exists, contents will be overwritten.")


    def getInfo(self):
        r = requests.get(self.recordID['url']['get'])

        if r.ok:
            print(f"\n*** Found Zenodo record {self.recordID['zenID']}: {r.json()['metadata']['title']}")
            print(f"Zenodo URL: {self.recordID['url']['doi']}")
            self.r = r

            self.downloadSize = sum(item['size'] for item in r.json()['files'])
            print(f"Record {self.recordID['zenID']}: {len(r.json()['files'])} files, {convert_bytes(self.downloadSize)}")

            # Display HTML job info from Zenodo page
            # This will render correctly in a notebook
            jobInfoStr = r.json()['metadata']['description']
            jobInfoHTML = HTML(jobInfoStr)
            display(jobInfoHTML)

            # Citation info link using hrefParser()
            p = hrefParser()
            p.feed(jobInfoStr)
            citeStr = f"Citation details: {p.output_list[0] + '#Cite-this-dataset'}"
            print(citeStr)
            self.recordID['url']['epsdata'] = p.output_list[0]

            # Set in object
            self.record = {'title':r.json()['metadata']['title'],
                           'jobInfo':jobInfoStr,
                           'jobInfoHTML':jobInfoHTML,
                           'cite':citeStr,
                           'ID':self.recordID}


        else:
            print(f"Can't find Zenodo record {recordID['zenID']}, error code: {r}")
            # self.r = r



    def downloadFiles(self, downloadList=[], overwriteFlag = False, overwritePromptFlag = True, hashCheck = False):
        """
        Download files from Zenodo record.

        Parameters
        ----------
        downloadList : list, optional, default=[]
            Set files to download from self.r.json()['files'] list.
            Defaults to all files.

        overwriteFlag : bool, optional, default = False
        overwritePromptFlag : bool, optional, default = True
            Options for reacquiring downloads if files already exist locally.

        TODO:

        - add hash check.

        """

        if not downloadList:
            downloadList = self.r.json()['files']

        fList=[]

        for n, item in enumerate(downloadList):
            print(f"\n***Getting item {item['links']['self']}")

            # Check if file exists locally already, and test size
            localFile = self.recordID['downloadDir']/item['key']

            downloadFlag = True

            # Slightly ugly logic here, but left all cases in for now for clarity while testing.
            if localFile.is_file():
                sizeCheck = localFile.stat().st_size - item['size']  # Quick file size check

                if sizeCheck:
                    print(f'Local file size incomensurate with remote by {sizeCheck} bytes. File will be downloaded again.')
                    downloadFlag = True

                else:
                    print('Local file already exists, file size OK.')

                    # if not (overwriteFlag and overwritePromptFlag): # Not correct - ends up as a catch-all apart from True-True case.
                    if not overwriteFlag:
                        downloadFlag = False

                    # elif (overwriteFlag and overwritePromptFlag):  # If true & true
                    else:
                        downloadFlag = True # Default case.

                        if overwritePromptFlag:
                            test = input("Download file again (y/n)?: ")

                            if test == 'y':
                                downloadFlag = True
                            else:
                                downloadFlag = False


                    # else:
                    #     downloadFlag = True

                    if downloadFlag:
                        print('File will be downloaded again.')

                        # Remove existing file in this case, otherwise wget will increment file name for new download.
                        print('Removing existing version...')
                        localFile.unlink()

                    else:
                        print('Skipping download.')

            if downloadFlag:
                # fout = wget.download(item['links']['self'], out=self.recordID['downloadDir'].as_posix())
                fout = wget.download(item['links']['self'], out=localFile.as_posix())

                print(f"Pulled to file: {fout}")
                fList.append(Path(fout))  # Log local file list

            else:
                # In this case, just rebuild local file list
                print(f"Existing file OK: {localFile}")
                fList.append(localFile)

            if hashCheck:
                print("Hash check not yet implemented.")


        self.fList = fList


    def unzipFiles(self, confirmFlag = True, verbose = True):
        # To do
        # - Add checking for multipart zip. Not sure of best way - from file suffix?
        # - Could use Counter maybe, as per epsman._repo.fileListCheck()
        # - Or look for unique file names as *.z*?
        #

        self.zip = {}
        # Unzip if required
        # for n, item in enumerate(self.fList):
            # if item.suffix == '.zip':
                # with zipfile.ZipFile(item,"r") as zipObj:
                #     zipFiles = zipObj.namelist()
                #     zipObj.extractall(self.recordID['downloadDir'])
                #
                #     print(f"Unzipped file {item}")
                #     self.zip[item.name]={'path':self.recordID['downloadDir'],
                #                          'zipfile':item,
                #                          'files':zipFiles}
                #
                # Version with local function
                # self.zip[item.name] = unzipFile(item)

        # With multipart zip checks
        # (1) List .z* suffixes, check for multipart as .z01
        zipList = []
        mpZipList = []
        [zipList.append(item) for item in self.fList if item.suffix.startswith('.zip')]
        [mpZipList.append(item) for item in self.fList if item.suffix.startswith('.z01')]

        print(f"\n*** Found {len(zipList)} archive(s).")

        if len(mpZipList) > 0:
            print("\n*** Found multipart archives, these will be joined & unzipped before full file extraction.")

            [joinZip(item.with_suffix('.zip')) for item in mpZipList]
            self.zipMP = [unzipFile(item.with_name(item.stem + '_joined.zip'), unzipPath=self.recordID['downloadDir']) for item in mpZipList]


            # The joined files contain the full zip file, so also need unzipping! Ugh.
            # Keep this as a loop just in case there are multiple .zip files (although there should only be one)
            for n, zipItem in enumerate(self.zipMP):

                # Fix issue with relative paths in some datasets - move file in this case.
                # Assume that extraneous paths include 'pkg', correlated with original pkg dir at Zenodo upload.

                # 1st attempt, with relative paths.
                # mvList = []
                # if 'pkg' in zipItem['zipfile'].relative_to(zipItem['path']).parts:
                    # shutil.move((zipItem['path']/zipItem['files'][0]).as_posix(), ABCOdata.zipMP[0]['path'].as_posix())
                    # mvList = [shutil.move((zipItem['path']/fileItem).as_posix(), zipItem['path'].as_posix()) for fileItem in zipItem['files']]
                    # TO DO: remove extraneous dir tree here too, not sure of best & safe way to do this!

                # Redo as loop over files - above doesn't work since relative paths not correct.
                mvList = []
                for fileItem in zipItem['files']:
                    fileItem = Path(fileItem)
                    if 'pkg' in fileItem.parts:
                        mvList.append(shutil.move((zipItem['path']/fileItem).as_posix(), zipItem['path'].as_posix()))

                        # Alt: could just use Path.rename() here too,
                        # e.g. testOut = (ABCOdata.zipMP[0]['path']/ABCOdata.zipMP[0]['files'][0]).rename(ABCOdata.zipMP[0]['path']/Path(ABCOdata.zipMP[0]['files'][0]).name)

                        # TO DO: remove extraneous dir tree here too, not sure of best & safe way to do this!
                        # Path.rmdir() should be safe way to do this - removes only empty dirs (but not a full tree).
                        # Ah, os.removedirs() is safe - stops at first non-empty dir in tree.
                        pruneDirTree((zipItem['path']/fileItem).parent)

                    else:
                        mvList.append(fileItem)

                # Hack to allow for cases where full path is included in zip file.
                # For linux zip this can be discarded with '-j', but doesn't seem to be an easy way in Python?
                # TODO: this still works, but can be tidied up now that mvList handles extraneous paths in multipart zips.
                localPaths = [Path(self.recordID['downloadDir'], item) for item in mvList]
                self.zip = [unzipFile(item, unzipPath=self.recordID['downloadDir'], confirmFlag = confirmFlag, verbose = verbose) for item in localPaths]

        else:
            self.zip = [unzipFile(item, unzipPath=self.recordID['downloadDir'], confirmFlag = confirmFlag, verbose = verbose) for item in zipList]


        # Alternative with counter - bit verbose. More robust?
        # [zipList.append(item) for item in self.fList if item.suffix.statswith('.z')]
        # # (2) Count file suffixes
        # suffixList = [Path(item).suffix for item in zipList]
        # c = Counter(suffixList)
        #
        # # (3) Check for multiple files
        # print(f"Found {c['.zip']} archives.")
        #
        # if c['.z01'] > 0:
        #     print("Found multipart archives, these will be joined before unzip.")

        # Run file sort routine to check and log.
        self.sortRecordFiles()



    def sortRecordFiles(self, verbose = True):
        '''
        Basic sorting/parsing of downloaded record.

        Based on epsman._repo.fileListCheck()

        '''

        for arch in self.zip:
            fileListTest = [self.record['ID']['downloadDir']/item for item in arch['files']]

            # Count file suffixes
            suffixList = [Path(item).suffix for item in fileListTest]
            c = Counter(suffixList)

            # Check for dirs & log - seems a bit redundant?
            # dirList = []
            # fileList = []
            # for n, item in enumerate(suffixList):
            #     # if Path(item).is_dir():       # Checking with Path OK for local files, but not for a list from another machine.
            #     if item == '':
            #         dirList.append(fileListTest[n])
            #     else:
            #         fileList.append(fileListTest[n])

            # Set ePS output file(s)
            self.ePSout = [item for item in fileListTest if item.suffix == '.out']

            # Set subdirs
            # Note these methods look at the full filesystem, not just arch files as per above.
            # Method with glob, could also use os.scandir, os.walk.
            # glob.glob(path.as_posix() + '/**/', recursive=True)
            # With os.walk
            self.record['files'] = list(os.walk(self.record['ID']['downloadDir']))
            self.record['dirList'] = [item[0] for item in self.record['files']]


            # Print details
            if verbose:
                print('\n***Record summary')
                print(f"Record {self.record['ID']['doi']}, title: {self.record['title']}")
                # print(f"Notebook file: {self.nbDetails[key]['file']}")
                # print(*self.nbDetails[key]['pkgInfo'][2:], sep='\n')
                # print(f"Root dir: {self.nbDetails[key]['pkgRootDir']}")

                print(f"Base dir: {self.record['ID']['downloadDir']}")
                print(f"Found {len(self.record['dirList'])-1} directories:")
                # print(*self.record['dirList'], sep='\n')
                [print('\t' + Path(item).relative_to(self.record['ID']['downloadDir']).as_posix()) for item in self.record['dirList'][1:]]

                print(f"\nFound {len(fileListTest)} items, with file types:")
                # Throw out full dictionary here, may want to switch to formatted list.
                # Set small width to force vertical format
                pprint.pprint(c, width=50)
