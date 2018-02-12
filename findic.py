import os

import astropy.io.fits as fits
import numpy as np
from dataanalysis import hashtools

import ast
import ddosa
import dataanalysis.core as da


class ICIndexEntry(ddosa.DataAnalysis):
    ds = "UNDEFINED"
    idx_hash="UNDEFINED"
    hashe = None
    version_from_index=True

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

class FindICIndexEntry(ddosa.DataAnalysis):
    ds = None
    icversion = 1
    version_from_index=True

    input_scw = ddosa.ScWData
    input_ic = ddosa.ICRoot

    run_for_hashe = True

    def find_entry(self):
        if hasattr(self,'input_scw') and hasattr(self.input_scw,'get_t1_t2'):
            t1, t2 = self.input_scw.get_t1_t2()
            revid=self.input_scw.input_scwid.str()[:4]
        else:
            t1=5000
            t2=5000
            revid="None"

        idxfn = self.input_ic.icroot + "/idx/ic/" + self.ds + "-IDX.fits"
        print("idx:", idxfn)

        idx_hash=hashtools.shhash(open(idxfn).read())[:8]

        idx = fits.open(idxfn)[1].data

        print(t1,t2,zip(idx['VSTART'],idx['VSTOP']))

        m_on = (idx['VSTART'] < t1) & (idx['VSTOP'] > t2) & (idx['VERSION'] == self.icversion)
        print("found valid:", sum(m_on))
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
                ic_hashe = ast.literal_eval(open(version_fn).read())

                print("searching",ic_hashe,rev_hashe)
                if hashtools.find_object(ic_hashe,rev_hashe):
                    ic_hashe_norev=hashtools.hashe_replace_object(ic_hashe,rev_hashe,'')
                    print("norev hashe",ic_hashe_norev)
                else:
                    print("unable to convert version file")

                return dict(hashe=ic_hashe, ds=self.ds, member_location=member_location, idx_hash=idx_hash, idxfn=idxfn)
            except IOError as e:
                print("unable to read version file, skipping",e)
            except SyntaxError as e:
                print("unable to read version file, skipping",e)

        print("found no version")
        return dict(hashe="UNDEFINED", ds=self.ds, member_location=member_location, idx_hash=idx_hash, idxfn=idxfn)

    def main(self):
        entry=self.find_entry()
        return ICIndexEntry(use_hashe=entry['hashe'], use_ds=entry['ds'], use_member_location=entry['member_location'], use_idx_hash=entry['idx_hash'], use_version_from_index=self.version_from_index)


