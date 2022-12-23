import os
import yaml

import astropy.io.fits as fits
import numpy as np
from dataanalysis import hashtools

import ast
import ddosa
import dataanalysis.core as da

class NoValidIC(ddosa.da.AnalysisException):
    pass

class ICIndexEntry(ddosa.DataAnalysis):
    ds = "UNDEFINED"
    idx_hash="UNDEFINED"
    hashe = None
    version_from_index=False

    def get_version(self):
        if self.version_from_index:
            v = self.get_signature() + "." + self.version + "." + self.ds + "-idx-" + self.idx_hash
        else:
            v = self.get_signature() + "." + self.version + "." + self.ds + "." + da.hashtools.shhash(self.hashe)[:8]

        generalized_hashe=self.hashe

        note = dict(
            origin_object=self.__class__.__name__,
            origin_module=__name__,
            generalized_hash=generalized_hashe,
            reduced_hash=v,
            handle=v,
        )

        if note not in self.factory.factorizations:
            self.factory.note_factorization(note)

        return v
        
    def get_member_location(self, scw=None):
        return self.member_location

       # icroot = self.icroot

        #icroot=os.environ['CURRENT_IC']
       # ic_version_code=icroot+"/ic_version.yaml"

       # if os.path.exists(ic_version_code):
       #     ic_version=yaml.load(open(ic_version_code))
      #  else:
       #     ic_version=dict(version_id=os.path.basename(icroot))

       # return ic_version['version_id']

class FindICIndexEntry(ddosa.DataAnalysis):
    ds = None
    icversion = 1
    version_from_index=True
    input_icroot=ddosa.ICRoot

 #   run_for_hashe = True
    #version_from_index=False

    def get_icversion(self, icroot):
        if "t20221103_osa11.2" not in icroot.ic_root_version:
            icversion = self.icversion
            print("using pre-set icversion", icversion)
        else:
            print("will use computed icversion")
            if (icversion:=getattr(self, '_icversion', None)) is not None:
                print("already computed icversion:", icversion)
            else:
                print("will derive icversion from memonic")
                master_file = fits.open(self.input_icroot.icroot + "/idx/ic/ic_master_file.fits")
                print(master_file[2].data)
                print(master_file[3].data)
                m_row = master_file[3].data[master_file[3].data['MNEMONIC'] == "OSA"]
                print("m_row", m_row)
                print("m_row", master_file[3].data.columns)
                self._icversion = m_row[self.ds.replace(".", "").replace("-", "_")]
                icversion = self._icversion

        return icversion


    def get_member_location(self,scw=None, icroot_obj=None):
        entry=self.find_entry(scw, icroot_obj=icroot_obj)
        return entry['member_location']

    def get_version(self):
        return self.get_signature()+"."+self.version #+"."+self.ic_version

    #@property
    #def icroot(self):
        #icroot=os.environ['CURRENT_IC']
        #i = ddosa.ICRoot()
        #i.main()
    #    return self.input_icroot.icroot

    def find_entry(self, scw, icroot_obj=None):

        if scw is not None:
            t1, t2 = scw.get_t1_t2()
            revid=scw.input_scwid.str()[:4]
        else:
            t1,t2=5000,5000
            revid="0000"

        print("calling find_entry with", scw, icroot_obj, self.input_icroot)
        icroot = (icroot_obj or self.input_icroot).icroot

        idxfn = icroot + "/idx/ic/" + self.ds + "-IDX.fits"
        print("idx:", idxfn)

        try:
            idx_hash=hashtools.shhash(open(idxfn, "rb").read())[:8]
        except Exception as e:
            print("problem reading this:", idxfn)
            raise

        idx = fits.open(idxfn)[1].data

        for v1, v2, vv in zip(idx['VSTART'],idx['VSTOP'],idx['VERSION']):
            print("requested:", t1, t2, self.get_icversion(icroot), "valid in IC:", v1, v2, vv)

        m_on = (idx['VSTART'] < t1) & (idx['VSTOP'] > t2) & (idx['VERSION'] == self.get_icversion(icroot))
        print("found valid:", sum(m_on))

        if sum(m_on)==0:
            raise NoValidIC("for range: %.10lg - %.10lg; have from %.10lg to %.10lg"%(t1, t2, min(idx['VSTART']), max(idx['VSTOP'])))

        print(idx[m_on])
        order = np.argsort(idx[m_on]['VSTART'])
        member_location_rel = idx[m_on]['MEMBER_LOCATION'][order[-1]]
        print("newest", member_location_rel)

        member_location = os.path.abspath(os.path.dirname(idxfn) + "/" + member_location_rel)

        assert os.path.exists(member_location)

        version_fn = os.path.dirname(member_location) + "/.version." + os.path.basename(member_location)

        print("version file", version_fn)

        rev_hashe=ddosa.Revolution(input_revid=revid).get_hashe()

        if os.path.exists(version_fn):
            print("found hashe file at", version_fn)

            try:
                version_f_content = open(version_fn).read()
                ic_hashe = ast.literal_eval(version_f_content)

                print("searching",ic_hashe,rev_hashe)
                if hashtools.find_object(ic_hashe,rev_hashe):
                    ic_hashe_norev=hashtools.hashe_replace_object(ic_hashe,rev_hashe,'')
                    print("norev hashe",ic_hashe_norev)
                else:
                    print("unable to extract revhashe ",rev_hashe,"from version file",ic_hashe)

                return dict(
                    hashe=ic_hashe, 
                    ds=self.ds, 
                    member_location=member_location, 
                    idx_hash=idx_hash, 
                    idxfn=idxfn,
                    icshhash=hashtools.shhash(ic_hashe)
                )
            except IOError as e:
                print("unable to read version file, skipping",e)
            except SyntaxError as e:
                print("unable to read version file, skipping",e)

        print("found no version")
        return dict(hashe="UNDEFINED", ds=self.ds, member_location=member_location, idx_hash=idx_hash, idxfn=idxfn)

    def main(self):
        pass
        #return ICIndexEntry(use_hashe=entry['hashe'], use_ds=entry['ds'], use_member_location=entry['member_location'], use_idx_hash=entry['idx_hash'], use_version_from_index=self.version_from_index)


