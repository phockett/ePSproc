"""
ePSproc ePSdata downloader

Basic IO function for downloading data from ePSdata, via Zenodo.

20/07/20    v1, basics in place.

TODO:

* More integration with ePSman. Use of json info files?
* Download multiple jobs/files. Could be in class, or as object per Zenodo record. Will be useful for, e.g., cases comparing multiple orbitals.
* Tidy up.

"""

#**** Imports
# Core
import requests
import wget
import zipfile
import os
from pathlib import Path

# Web stuff
from urllib.parse import urlparse
from IPython.core.display import HTML  # For HTML rendering support
from html.parser import HTMLParser

#**** Local functions
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
            print(f"/n*** Found Zenodo record {self.recordID['zenID']}: {r.json()['metadata']['title']}")
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
            print(f"Citation details: {p.output_list[0] + '#Cite-this-dataset'}")
            self.recordID['url']['epsdata'] = p.output_list[0]


        else:
            print(f"Can't find Zenodo record {recordID['zenID']}, error code: {r}")
            # self.r = r



    def downloadFiles(self, downloadList=[]):
        """
        Download files from Zenodo record.

        Parameters
        ----------
        downloadList : list, optional, default=[]
            Set files to download from self.r.json()['files'] list.
            Defaults to all files.

        """

        if not downloadList:
            downloadList = self.r.json()['files']

        fList=[]

        for n, item in enumerate(downloadList):
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
