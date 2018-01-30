import os

import astropy.io.fits as fits
import numpy as np

import ast
import ddosa
import dataanalysis.core as da


class ICIndexEntry(ddosa.DataAnalysis):
    ds = "UNDEFINED"
    hashe = None

    def get_version(self):
        v = self.get_signature() + "." + self.version + "." + self.ds + "." + da.hashtools.shhash(self.hashe)[:8]

        note = dict(
            origin_object=self.__class__.__name__,
            origin_module=__name__,
            generalized_hash=self.hashe,
            reduced_hash=v,
            handle=v,
        )

        if note not in self.factory.factorizations:
            self.factory.note_factorization(note)

        return v

class FindICIndexEntry(ddosa.DataAnalysis):
    ds = None
    icversion = 1
    input_scw = ddosa.ScWData
    input_ic = ddosa.ICRoot

    run_for_hashe = True

    def main(self):
        if hasattr(self,'input_scw'):
            t1, t2 = self.input_scw.get_t1_t2()
        else:
            t1=5000
            t2=5000

        idxfn = self.input_ic.icroot + "/idx/ic/" + self.ds + "-IDX.fits"
        print("idx:", idxfn)

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

        if os.path.exists(version_fn):
            print("found hashe file at", version_fn)

            try:
                ic_hashe = ast.literal_eval(open(version_fn).read())

                return ICIndexEntry(use_hashe=ic_hashe, use_ds=self.ds, use_member_location=member_location)
            except Exception as e:
                print("unable to read version file, skipping",e)

        return ICIndexEntry(use_hashe="UNDEFINED", use_ds=self.ds, use_member_location=member_location)
        #else:
        #    raise Exception("unable for find entry "+repr(self.ds)+" for "+repr(self.input_scw))
            # return DataAnalysis.from_hashe(ic_hashe).get()

