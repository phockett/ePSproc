"""
ePSproc ePSdata downloader

Basic IO function for downloading data from ePSdata, via Zenodo.

20/07/20    v1

"""

# Imports
import requests
import wget
import zipfile
import os
from pathlib import Path
from urllib.parse import urlparse

# Basic bytes to KB/Mb... conversion, from https://stackoverflow.com/questions/2104080/how-to-check-file-size-in-python
def convert_bytes(num):
    """
    This function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
#             return [num, x]
        num /= 1024.0


class ePSdata():
    """
    Class for interfacing with ePSdata via Zenodo API.

    Init by passing either DOI, URL or Zenodo ID.

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


    def setRecord(self, doi=None, ID=None, URL=None):
        # Set base locations - may want to do this elsewhere, and/or with looping dicts?
        doiBase = '10.5281/zenodo.'  # Common Zenodo DOI base
        doiURLBase = 'http://dx.doi.org/'
        zenURLBase = 'https://zenodo.org/record/'  # zenURLBase + ID corresponds to record webpage.
        zenAPIBase = 'https://zenodo.org/api/records/'  # Need /api for API calls!

        self.base = {}  # TODO - set above to dict item.

        # Set record IDs, starting from DOI
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
                recordID['doi'] = doiBase + ID
                # recordID['url']['doi'] = urljoin(doiURLbase + doiBase, ID)  # With urljoin
                recordID['url']['doi'] = doiURLbase + doiBase + ID  # From strs.


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
            self.recordID['downloadBase'] = Path(os.getcwd())
        else:
            self.recordID['downloadBase'] = Path(downloadBase)

        self.recordID['downloadDir'] = self.recordID['downloadBase']/str(self.recordID['zenID'])

        print(f"*** Download dir set to: {self.recordID['downloadDir']}")

    def makeLocal(self):
        try:
            os.mkdir(self.recordID['downloadDir'])
            print(f"Created {self.recordID['downloadDir']}")
        except FileExistsError:
            print(f"*** Directory {self.recordID['downloadDir']} already exists, contents will be overwritten.")

    def getInfo(self):
        r = requests.get(self.recordID['url']['get'])

        if r.ok:
            print(f"Found Zenodo record {self.recordID['zenID']}: {r.json()['metadata']['title']}")
            self.r = r

            self.downloadSize = sum(item['size'] for item in r.json()['files'])
            print(f"Record {self.recordID['zenID']}: {len(r.json()['files'])} files, {convert_bytes(self.downloadSize)}")

        else:
            print(f"Can't find Zenodo record {recordID['zenID']}, error code: {r}")
            # self.r = r



    def downloadFiles(self):

        fList = []

        for n, item in enumerate(self.r.json()['files']):
            print(f"Getting item {item['links']['self']}")
            fout = wget.download(item['links']['self'], out=self.recordID['downloadDir'].as_posix())

            print(f"Pulled to file: {fout}")
            fList.append(Path(fout))  # Log local file list

        self.fList = fList


    def unzipFiles(self):
        self.zip = {}
        # Unzip if required
        for n, item in enumerate(self.fList):
            if item.suffix == '.zip':
                with zipfile.ZipFile(item,"r") as zipObj:
                    zipFiles = zipObj.namelist()
                    zipObj.extractall(self.recordID['downloadDir'])

                    print(f"Unzipped file {item}")
                    self.zip[item.name]={'path':self.recordID['downloadDir'],
                                         'zipfile':item, 
                                         'files':zipFiles}
