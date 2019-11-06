#include <stdio.h>
extern "C" {
   int ncount_hepmc2=0 ;
   int hepmc2_delete_writer_(const int & position)
    {
        return 0;
        
    }
    int hepmc2_convert_event_(const int & position)
    {
        return 0;
    }
    int hepmc2_write_event_(const int & position)
    {
        if ( ncount_hepmc2 < 10) {
            printf(" RAPGAP: dummy version of write_event (HEPMC3) is used\n" );
            ++ncount_hepmc2;
        }
        return 0;
    }
    int hepmc2_clear_event_(const int & position)
    {
        return 0;
    }
    int hepmc2_set_cross_section_(const int & position, const double& x,const double& xe, const int& n1,const int& n2)
    {
        return 0;
    }
    int hepmc2_set_hepevt_address_(int* a)
    {
        return 0;
    }
    int hepmc2_set_attribute_int_(const int & position,const int & attval,const char* attname, size_t len)
    {
        return 0;
    }
    int hepmc2_set_attribute_double_(const int & position,const double & attval,const char* attname, size_t len)
    {
        return 0;
    }
    int hepmc2_new_writer_(const int & position,const int & mode,const char* ffilename, size_t len)
    {
        return 0;
    }
    int hepmc2_set_pdf_info_(const int & position,const int& parton_id1, const int& parton_id2, const double& x1, const double& x2,
                      const double& scale_in, const double& xf1,const double& xf2,
                      const int& pdf_id1, const int& pdf_id2)
    {
        return 0;
    }
    
}
