import numpy as np
from astropy.io import fits

def git_version():
    from subprocess import Popen, PIPE
    gitproc = Popen(['git', 'rev-parse','HEAD'], stdout = PIPE)
    (stdout, _) = gitproc.communicate()
    return stdout.strip()

def create_SPIRE_cat_from_Table(dataTable, prior250):
    """creates the XIDp catalogue in fits format required by HeDaM"""
    import datetime
    from astropy.table import Table

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='ID', format='100A', array=dataTable['HELP_ID'])
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=dataTable['RA'])
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=dataTable['Dec'])
    c4 = fits.Column(name='F_SPIRE_250', format='E', unit='mJy', array=dataTable['F_SPIRE_250'])
    c5 = fits.Column(name='FErr_SPIRE_250_u', format='E', unit='mJy', array=dataTable['FErr_SPIRE_250_u'])
    c6 = fits.Column(name='FErr_SPIRE_250_l', format='E', unit='mJy', array=dataTable['FErr_SPIRE_250_l'])
    c7 = fits.Column(name='F_SPIRE_350', format='E', unit='mJy', array=dataTable['F_SPIRE_350'])
    c8 = fits.Column(name='FErr_SPIRE_350_u', format='E', unit='mJy', array=dataTable['FErr_SPIRE_350_u'])
    c9 = fits.Column(name='FErr_SPIRE_350_l', format='E', unit='mJy', array=dataTable['FErr_SPIRE_350_l'])
    c10 = fits.Column(name='F_SPIRE_500', format='E', unit='mJy', array=dataTable['F_SPIRE_500'])
    c11 = fits.Column(name='FErr_SPIRE_500_u', format='E', unit='mJy', array=dataTable['FErr_SPIRE_500_u'])
    c12 = fits.Column(name='FErr_SPIRE_500_l', format='E', unit='mJy', array=dataTable['FErr_SPIRE_500_l'])
    c13 = fits.Column(name='Bkg_SPIRE_250', format='E', unit='mJy', array=dataTable['Bkg_SPIRE_250'])
    c14 = fits.Column(name='Bkg_SPIRE_350', format='E', unit='mJy', array=dataTable['Bkg_SPIRE_350'])
    c15 = fits.Column(name='Bkg_SPIRE_500', format='E', unit='mJy', array=dataTable['Bkg_SPIRE_500'])
    c16 = fits.Column(name='Sig_conf_SPIRE_250', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_SPIRE_250'])
    c17 = fits.Column(name='Sig_conf_SPIRE_350', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_SPIRE_350'])
    c18 = fits.Column(name='Sig_conf_SPIRE_500', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_SPIRE_500'])
    c19 = fits.Column(name='Rhat_SPIRE_250', format='E', array=dataTable['Rhat_SPIRE_250'])
    c20 = fits.Column(name='Rhat_SPIRE_350', format='E', array=dataTable['Rhat_SPIRE_350'])
    c21 = fits.Column(name='Rhat_SPIRE_500', format='E', array=dataTable['Rhat_SPIRE_500'])
    c22 = fits.Column(name='n_eff_SPIRE_250', format='E', array=dataTable['n_eff_SPIRE_250'])
    c23 = fits.Column(name='n_eff_SPIRE_350', format='E', array=dataTable['n_eff_SPIRE_350'])
    c24 = fits.Column(name='n_eff_SPIRE_500', format='E', array=dataTable['n_eff_SPIRE_500'])
    c25 = fits.Column(name='Pval_res_250', format='E', array=dataTable['Pval_res_250'])
    c26 = fits.Column(name='Pval_res_350', format='E', array=dataTable['Pval_res_250'])
    c27 = fits.Column(name='Pval_res_500', format='E', array=dataTable['Pval_res_250'])
    c28 = fits.Column(name='tile', format='K', array=dataTable['tile'].astype(float).astype(int))

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,
                                           c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28])

    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2')

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3')

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')
    tbhdu.header.set('TDESC4','250 Flux (at 50th percentile)',after='TUCD4')

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')
    tbhdu.header.set('TDESC5','250 Flux (at 84.1 percentile) ',after='TUCD5')

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')
    tbhdu.header.set('TDESC6','250 Flux (at 15.9 percentile)',after='TUCD6')

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')
    tbhdu.header.set('TDESC7','350 Flux (at 50th percentile)',after='TUCD7')

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')
    tbhdu.header.set('TDESC8','350 Flux (at 84.1 percentile) ',after='TUCD8')

    tbhdu.header.set('TUCD9','phot.flux.density',after='TUNIT9')
    tbhdu.header.set('TDESC9','350 Flux (at 15.9 percentile)',after='TUCD9')

    tbhdu.header.set('TUCD10','phot.flux.density',after='TUNIT10')
    tbhdu.header.set('TDESC10','500 Flux (at 50th percentile)',after='TUCD10')

    tbhdu.header.set('TUCD11','phot.flux.density',after='TUNIT11')
    tbhdu.header.set('TDESC11','500 Flux (at 84.1 percentile) ',after='TUCD11')

    tbhdu.header.set('TUCD12','phot.flux.density',after='TUNIT12')
    tbhdu.header.set('TDESC12','500 Flux (at 15.9 percentile)',after='TUCD12')

    tbhdu.header.set('TUCD13','phot.flux.density',after='TUNIT13')
    tbhdu.header.set('TDESC13','250 background',after='TUCD13')

    tbhdu.header.set('TUCD14','phot.flux.density',after='TUNIT14')
    tbhdu.header.set('TDESC14','350 background',after='TUCD14')

    tbhdu.header.set('TUCD15','phot.flux.density',after='TUNIT15')
    tbhdu.header.set('TDESC15','500 background',after='TUCD15')

    tbhdu.header.set('TUCD16','phot.flux.density',after='TUNIT16')
    tbhdu.header.set('TDESC16','250 residual confusion noise',after='TUCD16')

    tbhdu.header.set('TUCD17','phot.flux.density',after='TUNIT17')
    tbhdu.header.set('TDESC17','350 residual confusion noise',after='TUCD17')

    tbhdu.header.set('TUCD18','phot.flux.density',after='TUNIT18')
    tbhdu.header.set('TDESC18','500 residual confusion noise',after='TUCD18')

    tbhdu.header.set('TUCD19','stat.value',after='TFORM16')
    tbhdu.header.set('TDESC19','250 MCMC Convergence statistic',after='TUCD19')

    tbhdu.header.set('TUCD20','stat.value',after='TFORM20')
    tbhdu.header.set('TDESC20','350 MCMC Convergence statistic',after='TUCD20')

    tbhdu.header.set('TUCD21','stat.value',after='TFORM21')
    tbhdu.header.set('TDESC21','500 MCMC Convergence statistic',after='TUCD21')

    tbhdu.header.set('TUCD22','stat.value',after='TFORM22')
    tbhdu.header.set('TDESC22','250 MCMC independence statistic',after='TUCD22')

    tbhdu.header.set('TUCD23','stat.value',after='TFORM23')
    tbhdu.header.set('TDESC23','350 MCMC independence statistic',after='TUCD23')

    tbhdu.header.set('TUCD24','stat.value',after='TFORM24')
    tbhdu.header.set('TDESC24','500 MCMC independence statistic',after='TUCD24')
    
    tbhdu.header.set('TUCD25','stat.value',after='TFORM25')
    tbhdu.header.set('TDESC25','250 Bayes Pval residual statistic',after='TUCD25')

    tbhdu.header.set('TUCD26','stat.value',after='TFORM26')
    tbhdu.header.set('TDESC26','350 Bayes Pval residual statistic',after='TUCD26')

    tbhdu.header.set('TUCD27','stat.value',after='TFORM27')
    tbhdu.header.set('TDESC27','500 Bayes Pval residual statistic',after='TUCD27')

    tbhdu.header.set('TUCD28','value',after='TFORM28')
    tbhdu.header.set('TDESC28','Tile Number',after='TUCD28')

    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Pri_Cat'] = prior250.prior_cat_file
    prihdr['TITLE']   = 'SPIRE XID+ catalogue'
    #prihdr['OBJECT']  = prior250.imphdu['OBJECT'] #I need to think if this needs to change
    prihdr['CREATOR'] = 'WP5'
    #prihdr['XID+V'] = git_version()
    prihdr['DATE']    = datetime.datetime.now().isoformat()
    prihdu = fits.PrimaryHDU(header=prihdr)

    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist

def make_master_SPIRE_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=0):
    """Make master source catalogue from XID+ with HEALpix
       Master_filename needs the .pkl extension
       Creates Master_Catalogue.fits in output_folder
       Will overwrite existing Master_Catalogue.fits if one exists"""

    import pickle
    import dill
    from astropy.table import Table, vstack
    from xidplus import moc_routines, catalogue

    with open(output_folder+Master_filename, "rb") as f:
        obj=pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']
    prior250 = obj['psw']

    #Need a master catalogue to add to so do the 0th tile outside the loop
    for i in range(start_tile, len(tiles)):
        print('On tile '+str(i)+' of ' + str(len(tiles)))

        tile = tiles[i]
        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing?')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data?')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()
                main_start = i+1
                break

    tilepsw = obj['psw']
    tilepmw = obj['pmw']
    tileplw = obj['plw']
    tilepos = obj['posterior']

    master_cat = catalogue.create_SPIRE_cat(tilepos, tilepsw, tilepmw, tileplw)
    #Only keep sources inside the tile
    kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
    kept_sources = np.array(kept_sources)

    master_cat[1].data = master_cat[1].data[kept_sources]
    #Create the master catalogue table
    dataTable = Table(master_cat[1].data)
    dataTable['tile'] = main_start-1

    for i in range(main_start, len(tiles)):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)))
        tile = tiles[i]

        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()

                tilepsw = obj['psw']
                tilepmw = obj['pmw']
                tileplw = obj['plw']
                tilepos = obj['posterior']

                hdulist = catalogue.create_SPIRE_cat(tilepos, tilepsw, tilepmw, tileplw)
                #Only keep sources inside the tile
                kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
                kept_sources = np.array(kept_sources)

                hdulist[1].data = hdulist[1].data[kept_sources]
                tileTable = Table(hdulist[1].data)
                tileTable['tile'] = np.full(len(tileTable),i)
                #Append the tile catalogue to the master catalogue
                dataTable = vstack([dataTable, tileTable])

    #Write the master catalogue to fits
    master_cat = create_SPIRE_cat_from_Table(dataTable, prior250)
    master_cat.writeto(output_folder+'Master_Catalogue_SPIRE.fits', overwrite = True)

def make_master_SPIRE_catalogue_HEALpix_split(output_folder, Master_filename, tile_file_name, start_tile=0, end_tile=1e4):
    """Make master source catalogue from XID+ with HEALpix
       Master_filename needs the .pkl extension
       Creates Master_Catalogue.fits in output_folder
       Will overwrite existing Master_Catalogue.fits if one exists"""

    import pickle
    import dill
    from astropy.table import Table, vstack
    from xidplus import moc_routines, catalogue

    with open(output_folder+Master_filename, "rb") as f:
        obj=pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']
    prior250 = obj['psw']

    if end_tile > len(tiles):
        end_tile = len(tiles)

    #Need a master catalogue to add to so do the 0th tile outside the loop
    for i in range(start_tile, end_tile):
        print('On tile '+str(i)+' of ' + str(len(tiles)))

        tile = tiles[i]
        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing?')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data?')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()
                main_start = i+1
                break

    tilepsw = obj['psw']
    tilepmw = obj['pmw']
    tileplw = obj['plw']
    tilepos = obj['posterior']

    master_cat = catalogue.create_SPIRE_cat(tilepos, tilepsw, tilepmw, tileplw)
    #Only keep sources inside the tile
    kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
    kept_sources = np.array(kept_sources)

    master_cat[1].data = master_cat[1].data[kept_sources]
    #Create the master catalogue table
    dataTable = Table(master_cat[1].data)
    dataTable['tile'] = main_start-1

    for i in range(main_start, end_tile):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)))
        tile = tiles[i]

        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()

                tilepsw = obj['psw']
                tilepmw = obj['pmw']
                tileplw = obj['plw']
                tilepos = obj['posterior']

                hdulist = catalogue.create_SPIRE_cat(tilepos, tilepsw, tilepmw, tileplw)
                #Only keep sources inside the tile
                kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
                kept_sources = np.array(kept_sources)

                hdulist[1].data = hdulist[1].data[kept_sources]
                tileTable = Table(hdulist[1].data)
                tileTable['tile'] = np.full(len(tileTable),i)
                #Append the tile catalogue to the master catalogue
                dataTable = vstack([dataTable, tileTable])

    #Write the master catalogue to fits
    master_cat = create_SPIRE_cat_from_Table(dataTable, prior250)
    master_cat.writeto(output_folder+'Master_Catalogue_SPIRE_'+str(start_tile)+'.fits', overwrite = True)
    
def create_SPIRE_cat_post(posterior,prior250,prior350,prior500):
    """
    Create SPIRE catalogue of posterior

    :param posterior: SPIRE xidplus.posterior class
    :param prior250: SPIRE 250 xidplus.prior class
    :param prior350: SPIRE 350 xidplus.prior class
    :param prior500: SPIRE 500 xidplus.prior class
    :return: fits hdulist
    """
    import datetime
    nsrc=posterior.nsrc
    samp_chains=len(posterior.samples['src_f'])

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='ID', format='100A', array=prior250.ID)
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=prior250.sra)
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=prior250.sdec)
    c4 = fits.Column(name='F_SPIRE_250', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],50.0,axis=0))
    c5 = fits.Column(name='FErr_SPIRE_250_u', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],84.1,axis=0))
    c6 = fits.Column(name='FErr_SPIRE_250_l', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],15.9,axis=0))
    c7 = fits.Column(name='F_SPIRE_350', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,1,:],50.0,axis=0))
    c8 = fits.Column(name='FErr_SPIRE_350_u', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,1,:],84.1,axis=0))
    c9 = fits.Column(name='FErr_SPIRE_350_l', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,1,:],15.9,axis=0))
    c10 = fits.Column(name='F_SPIRE_500', format='E', unit='mJy',
                      array=np.percentile(posterior.samples['src_f'][:,2,:],50.0,axis=0))
    c11 = fits.Column(name='FErr_SPIRE_500_u', format='E', unit='mJy',
                      array=np.percentile(posterior.samples['src_f'][:,2,:],84.1,axis=0))
    c12 = fits.Column(name='FErr_SPIRE_500_l', format='E', unit='mJy',
                      array=np.percentile(posterior.samples['src_f'][:,2,:],15.9,axis=0))
    c13 = fits.Column(name='Bkg_SPIRE_250', format='E', unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['bkg'][:,0],50.0,axis=0)))
    c14 = fits.Column(name='Bkg_SPIRE_350', format='E', unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['bkg'][:,1],50.0,axis=0)))
    c15 = fits.Column(name='Bkg_SPIRE_500', format='E', unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['bkg'][:,2],50.0,axis=0)))
    c16 = fits.Column(name='Sig_conf_SPIRE_250', format='E',unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['sigma_conf'][:,0],50.0,axis=0)))
    c17 = fits.Column(name='Sig_conf_SPIRE_350', format='E',unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['sigma_conf'][:,1],50.0,axis=0)))
    c18 = fits.Column(name='Sig_conf_SPIRE_500', format='E',unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['sigma_conf'][:,2],50.0,axis=0)))
    c19 = fits.Column(name='Rhat_SPIRE_250', format='E', array=posterior.Rhat['src_f'][:,0])
    c20 = fits.Column(name='Rhat_SPIRE_350', format='E', array=posterior.Rhat['src_f'][:,1])
    c21 = fits.Column(name='Rhat_SPIRE_500', format='E', array=posterior.Rhat['src_f'][:,2])
    c22 = fits.Column(name='n_eff_SPIRE_250', format='E', array=posterior.n_eff['src_f'][:,0])
    c23 = fits.Column(name='n_eff_SPIRE_350', format='E', array=posterior.n_eff['src_f'][:,1])
    c24 = fits.Column(name='n_eff_SPIRE_500', format='E', array=posterior.n_eff['src_f'][:,2])
    c25 = fits.Column(name='Post_SPIRE_250', format=str(samp_chains)+'E', array=posterior.samples['src_f'][:,0,:].T)
    c26 = fits.Column(name='Post_SPIRE_350', format=str(samp_chains)+'E', array=posterior.samples['src_f'][:,1,:].T)
    c27 = fits.Column(name='Post_SPIRE_500', format=str(samp_chains)+'E', array=posterior.samples['src_f'][:,2,:].T)
    #c28 = fits.Column(name='Post_bkg_SPIRE_250', format=str(samp_chains)+'E', unit='mJy', array=np.repeat(posterior.samples['bkg'][:,0], nsrc, axis=0))
    #c29 = fits.Column(name='Post_bkg_SPIRE_350', format=str(samp_chains)+'E', unit='mJy', array=np.repeat(posterior.samples['bkg'][:,1], nsrc, axis=0))
    #c30 = fits.Column(name='Post_bkg_SPIRE_500', format=str(samp_chains)+'E', unit='mJy', array=np.repeat(posterior.samples['bkg'][:,2], nsrc, axis=0))

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27]) #,c28,c29,c30])
    
    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')      
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2') 

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')      
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3') 

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')      
    tbhdu.header.set('TDESC4','250 Flux (at 50th percentile)',after='TUCD4') 

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')      
    tbhdu.header.set('TDESC5','250 Flux (at 84.1 percentile) ',after='TUCD5') 

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')      
    tbhdu.header.set('TDESC6','250 Flux (at 15.9 percentile)',after='TUCD6') 

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')      
    tbhdu.header.set('TDESC7','350 Flux (at 50th percentile)',after='TUCD7') 

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')      
    tbhdu.header.set('TDESC8','350 Flux (at 84.1 percentile) ',after='TUCD8') 

    tbhdu.header.set('TUCD9','phot.flux.density',after='TUNIT9')      
    tbhdu.header.set('TDESC9','350 Flux (at 15.9 percentile)',after='TUCD9') 

    tbhdu.header.set('TUCD10','phot.flux.density',after='TUNIT10')      
    tbhdu.header.set('TDESC10','500 Flux (at 50th percentile)',after='TUCD10') 

    tbhdu.header.set('TUCD11','phot.flux.density',after='TUNIT11')      
    tbhdu.header.set('TDESC11','500 Flux (at 84.1 percentile) ',after='TUCD11') 

    tbhdu.header.set('TUCD12','phot.flux.density',after='TUNIT12')      
    tbhdu.header.set('TDESC12','500 Flux (at 15.9 percentile)',after='TUCD12')

    tbhdu.header.set('TUCD13','phot.flux.density',after='TUNIT13')      
    tbhdu.header.set('TDESC13','250 background',after='TUCD13') 

    tbhdu.header.set('TUCD14','phot.flux.density',after='TUNIT14')      
    tbhdu.header.set('TDESC14','350 background',after='TUCD14') 

    tbhdu.header.set('TUCD15','phot.flux.density',after='TUNIT15')      
    tbhdu.header.set('TDESC15','500 background',after='TUCD15')

    tbhdu.header.set('TUCD16','phot.flux.density',after='TUNIT16')
    tbhdu.header.set('TDESC16','250 residual confusion noise',after='TUCD16')

    tbhdu.header.set('TUCD17','phot.flux.density',after='TUNIT17')
    tbhdu.header.set('TDESC17','350 residual confusion noise',after='TUCD17')

    tbhdu.header.set('TUCD18','phot.flux.density',after='TUNIT18')
    tbhdu.header.set('TDESC18','500 residual confusion noise',after='TUCD18')

    tbhdu.header.set('TUCD19','stat.value',after='TFORM19')      
    tbhdu.header.set('TDESC19','250 MCMC Convergence statistic',after='TUCD19')

    tbhdu.header.set('TUCD20','stat.value',after='TFORM20')      
    tbhdu.header.set('TDESC20','350 MCMC Convergence statistic',after='TUCD20')

    tbhdu.header.set('TUCD21','stat.value',after='TFORM21')      
    tbhdu.header.set('TDESC21','500 MCMC Convergence statistic',after='TUCD21')

    tbhdu.header.set('TUCD22','stat.value',after='TFORM22')      
    tbhdu.header.set('TDESC22','250 MCMC independence statistic',after='TUCD22')

    tbhdu.header.set('TUCD23','stat.value',after='TFORM23')      
    tbhdu.header.set('TDESC23','350 MCMC independence statistic',after='TUCD23')

    tbhdu.header.set('TUCD24','stat.value',after='TFORM24')      
    tbhdu.header.set('TDESC24','500 MCMC independence statistic',after='TUCD24')

    tbhdu.header.set('TUCD25','phot.flux.density',after='TFORM25')      
    tbhdu.header.set('TDESC25','250 samples',after='TUCD25')

    tbhdu.header.set('TUCD26','phot.flux.density',after='TFORM26')      
    tbhdu.header.set('TDESC26','350 samples',after='TUCD26')
    
    tbhdu.header.set('TUCD27','phot.flux.density',after='TFORM27')      
    tbhdu.header.set('TDESC27','500 samples',after='TUCD27')

    #tbhdu.header.set('TUCD28','phot.flux.density',after='TFORM28')      
    #tbhdu.header.set('TDESC28','250 bkg samples',after='TUCD28')

    #tbhdu.header.set('TUCD29','phot.flux.density',after='TFORM29')      
    #tbhdu.header.set('TDESC29','350 bkg samples',after='TUCD29')
    
    #tbhdu.header.set('TUCD30','phot.flux.density',after='TFORM30')      
    #tbhdu.header.set('TDESC30','500 bkg samples',after='TUCD30')


    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Prior_C'] = prior250.prior_cat_file
    prihdr['TITLE']   = 'SPIRE XID catalogue'        
    #prihdr['OBJECT']  = prior250.imphdu['OBJECT'] #I need to think if this needs to change                              
    prihdr['CREATOR'] = 'WP5'                                 
    #prihdr['VERSION'] = 'beta'                                 
    prihdr['DATE']    = datetime.datetime.now().isoformat()              
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist

def create_SPIRE_cat_post_from_Table(dataTable, prior250, samples_chains):
    """creates the XIDp catalogue in fits format required by HeDaM"""
    import datetime
    from astropy.table import Table

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='ID', format='100A', array=dataTable['HELP_ID'])
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=dataTable['RA'])
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=dataTable['Dec'])
    c4 = fits.Column(name='F_SPIRE_250', format='E', unit='mJy', array=dataTable['F_SPIRE_250'])
    c5 = fits.Column(name='FErr_SPIRE_250_u', format='E', unit='mJy', array=dataTable['FErr_SPIRE_250_u'])
    c6 = fits.Column(name='FErr_SPIRE_250_l', format='E', unit='mJy', array=dataTable['FErr_SPIRE_250_l'])
    c7 = fits.Column(name='F_SPIRE_350', format='E', unit='mJy', array=dataTable['F_SPIRE_350'])
    c8 = fits.Column(name='FErr_SPIRE_350_u', format='E', unit='mJy', array=dataTable['FErr_SPIRE_350_u'])
    c9 = fits.Column(name='FErr_SPIRE_350_l', format='E', unit='mJy', array=dataTable['FErr_SPIRE_350_l'])
    c10 = fits.Column(name='F_SPIRE_500', format='E', unit='mJy', array=dataTable['F_SPIRE_500'])
    c11 = fits.Column(name='FErr_SPIRE_500_u', format='E', unit='mJy', array=dataTable['FErr_SPIRE_500_u'])
    c12 = fits.Column(name='FErr_SPIRE_500_l', format='E', unit='mJy', array=dataTable['FErr_SPIRE_500_l'])
    c13 = fits.Column(name='Bkg_SPIRE_250', format='E', unit='mJy', array=dataTable['Bkg_SPIRE_250'])
    c14 = fits.Column(name='Bkg_SPIRE_350', format='E', unit='mJy', array=dataTable['Bkg_SPIRE_350'])
    c15 = fits.Column(name='Bkg_SPIRE_500', format='E', unit='mJy', array=dataTable['Bkg_SPIRE_500'])
    c16 = fits.Column(name='Sig_conf_SPIRE_250', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_SPIRE_250'])
    c17 = fits.Column(name='Sig_conf_SPIRE_350', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_SPIRE_350'])
    c18 = fits.Column(name='Sig_conf_SPIRE_500', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_SPIRE_500'])
    c19 = fits.Column(name='Rhat_SPIRE_250', format='E', array=dataTable['Rhat_SPIRE_250'])
    c20 = fits.Column(name='Rhat_SPIRE_350', format='E', array=dataTable['Rhat_SPIRE_350'])
    c21 = fits.Column(name='Rhat_SPIRE_500', format='E', array=dataTable['Rhat_SPIRE_500'])
    c22 = fits.Column(name='n_eff_SPIRE_250', format='E', array=dataTable['n_eff_SPIRE_250'])
    c23 = fits.Column(name='n_eff_SPIRE_350', format='E', array=dataTable['n_eff_SPIRE_350'])
    c24 = fits.Column(name='n_eff_SPIRE_500', format='E', array=dataTable['n_eff_SPIRE_500'])
    c25 = fits.Column(name='Post_SPIRE_250', format=str(samples_chains)+'E', unit='mJy', array=dataTable['Post_SPIRE_250'])
    c26 = fits.Column(name='Post_SPIRE_350', format=str(samples_chains)+'E', unit='mJy', array=dataTable['Post_SPIRE_350'])
    c27 = fits.Column(name='Post_SPIRE_500', format=str(samples_chains)+'E', unit='mJy', array=dataTable['Post_SPIRE_500'])
    #c28 = fits.Column(name='Post_bkg_SPIRE_250', format=str(samples*chains)+'E', unit='mJy', array=dataTable['Post_bkg_SPIRE_250'])
    #c29 = fits.Column(name='Post_bkg_SPIRE_350', format=str(samples*chains)+'E', unit='mJy', array=dataTable['Post_bkg_SPIRE_350'])
    #c30 = fits.Column(name='Post_bkg_SPIRE_500', format=str(samples*chains)+'E', unit='mJy', array=dataTable['Post_bkg_SPIRE_500'])

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27]) #,c28,c29,c30])

    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2')

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3')

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')
    tbhdu.header.set('TDESC4','250 Flux (at 50th percentile)',after='TUCD4')

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')
    tbhdu.header.set('TDESC5','250 Flux (at 84.1 percentile) ',after='TUCD5')

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')
    tbhdu.header.set('TDESC6','250 Flux (at 15.9 percentile)',after='TUCD6')

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')
    tbhdu.header.set('TDESC7','350 Flux (at 50th percentile)',after='TUCD7')

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')
    tbhdu.header.set('TDESC8','350 Flux (at 84.1 percentile) ',after='TUCD8')

    tbhdu.header.set('TUCD9','phot.flux.density',after='TUNIT9')
    tbhdu.header.set('TDESC9','350 Flux (at 15.9 percentile)',after='TUCD9')

    tbhdu.header.set('TUCD10','phot.flux.density',after='TUNIT10')
    tbhdu.header.set('TDESC10','500 Flux (at 50th percentile)',after='TUCD10')

    tbhdu.header.set('TUCD11','phot.flux.density',after='TUNIT11')
    tbhdu.header.set('TDESC11','500 Flux (at 84.1 percentile) ',after='TUCD11')

    tbhdu.header.set('TUCD12','phot.flux.density',after='TUNIT12')
    tbhdu.header.set('TDESC12','500 Flux (at 15.9 percentile)',after='TUCD12')

    tbhdu.header.set('TUCD13','phot.flux.density',after='TUNIT13')
    tbhdu.header.set('TDESC13','250 background',after='TUCD13')

    tbhdu.header.set('TUCD14','phot.flux.density',after='TUNIT14')
    tbhdu.header.set('TDESC14','350 background',after='TUCD14')

    tbhdu.header.set('TUCD15','phot.flux.density',after='TUNIT15')
    tbhdu.header.set('TDESC15','500 background',after='TUCD15')

    tbhdu.header.set('TUCD16','phot.flux.density',after='TUNIT16')
    tbhdu.header.set('TDESC16','250 residual confusion noise',after='TUCD16')

    tbhdu.header.set('TUCD17','phot.flux.density',after='TUNIT17')
    tbhdu.header.set('TDESC17','350 residual confusion noise',after='TUCD17')

    tbhdu.header.set('TUCD18','phot.flux.density',after='TUNIT18')
    tbhdu.header.set('TDESC18','500 residual confusion noise',after='TUCD18')

    tbhdu.header.set('TUCD19','stat.value',after='TFORM16')
    tbhdu.header.set('TDESC19','250 MCMC Convergence statistic',after='TUCD19')

    tbhdu.header.set('TUCD20','stat.value',after='TFORM20')
    tbhdu.header.set('TDESC20','350 MCMC Convergence statistic',after='TUCD20')

    tbhdu.header.set('TUCD21','stat.value',after='TFORM21')
    tbhdu.header.set('TDESC21','500 MCMC Convergence statistic',after='TUCD21')

    tbhdu.header.set('TUCD22','stat.value',after='TFORM22')
    tbhdu.header.set('TDESC22','250 MCMC independence statistic',after='TUCD22')

    tbhdu.header.set('TUCD23','stat.value',after='TFORM23')
    tbhdu.header.set('TDESC23','350 MCMC independence statistic',after='TUCD23')

    tbhdu.header.set('TUCD24','stat.value',after='TFORM24')
    tbhdu.header.set('TDESC24','500 MCMC independence statistic',after='TUCD24')

    tbhdu.header.set('TUCD25','phot.flux.density',after='TFORM25')
    tbhdu.header.set('TDESC25','250 samples',after='TUCD25')

    tbhdu.header.set('TUCD26','phot.flux.density',after='TFORM26')
    tbhdu.header.set('TDESC26','350 samples',after='TUCD26')

    tbhdu.header.set('TUCD27','phot.flux.density',after='TFORM27')
    tbhdu.header.set('TDESC27','500 samples',after='TUCD27')

    #tbhdu.header.set('TUCD28','phot.flux.density',after='TFORM28')
    #tbhdu.header.set('TDESC28','250 bkg samples',after='TUCD28')

    #tbhdu.header.set('TUCD29','phot.flux.density',after='TFORM29')
    #tbhdu.header.set('TDESC29','350 bkg samples',after='TUCD29')

    #tbhdu.header.set('TUCD30','phot.flux.density',after='TFORM30')
    #tbhdu.header.set('TDESC30','500 bkg samples',after='TUCD30')

    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Pri_Cat'] = prior250.prior_cat_file
    prihdr['TITLE']   = 'SPIRE XID+ Posterior catalogue'
    #prihdr['OBJECT']  = prior250.imphdu['OBJECT'] #I need to think if this needs to change
    prihdr['CREATOR'] = 'WP5'
    #prihdr['XID+V'] = git_version()
    prihdr['DATE']    = datetime.datetime.now().isoformat()
    prihdu = fits.PrimaryHDU(header=prihdr)

    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist

def make_master_SPIRE_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=0, end_tile=0):
    """Make master source catalogue from XID+ with HEALpix
       Master_filename needs the .pkl extension
       Creates Master_Catalogue.fits in output_folder
       Will overwrite existing Master_Catalogue.fits if one exists"""

    import pickle
    import dill
    from astropy.table import Table, vstack
    from xidplus import moc_routines, catalogue

    with open(output_folder+Master_filename, "rb") as f:
        obj=pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']
    prior250 = obj['psw']

    if end_tile == 0:
        end_tile = len(tiles)

    #Need a master catalogue to add to so do the 0th tile outside the loop
    for i in range(start_tile, end_tile):
        print('On tile '+str(i)+' of ' + str(len(tiles)))

        tile = tiles[i]
        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing?')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data?')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()
                main_start = i+1
                break

    tilepsw = obj['psw']
    tilepmw = obj['pmw']
    tileplw = obj['plw']
    tilepos = obj['posterior']

    master_cat = create_SPIRE_cat_post(tilepos, tilepsw, tilepmw, tileplw)
    #Only keep sources inside the tile
    kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
    kept_sources = np.array(kept_sources)

    master_cat[1].data = master_cat[1].data[kept_sources]
    #Create the master catalogue table
    dataTable = Table(master_cat[1].data)
    dataTable['tile'] = main_start-1

    for i in range(main_start, end_tile):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)))
        tile = tiles[i]

        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()

                tilepsw = obj['psw']
                tilepmw = obj['pmw']
                tileplw = obj['plw']
                tilepos = obj['posterior']

                hdulist = create_SPIRE_cat_post(tilepos, tilepsw, tilepmw, tileplw)
                #Only keep sources inside the tile
                kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
                kept_sources = np.array(kept_sources)

                hdulist[1].data = hdulist[1].data[kept_sources]
                tileTable = Table(hdulist[1].data)
                #Append the tile catalogue to the master catalogue
                dataTable = vstack([dataTable, tileTable])

    #Write the master catalogue to fits
    samples_chains = len(tilepos.samples['src_f'])
    master_cat = create_SPIRE_cat_post_from_Table(dataTable, prior250, samples_chains)
    master_cat.writeto(output_folder+'Master_Post_Catalogue_SPIRE_'+str(start_tile)+'.fits', overwrite = True)
    

def make_master_map_HEALpix(output_folder, Master_filename, tile_file_name, band, im, imhdu, w_pri, start_tile=0):
    """Create a replicated map from XID+ with HEALpix tiles
       Will overwrite existing Master_Map_###.fits if one exists"""

    import pickle
    from scipy.sparse import coo_matrix
    from astropy import wcs
    from xidplus import moc_routines as moc
    
    sam = 2999
    
    if band == 'psw':
        b = 0
    elif band == 'pmw':
        b = 1
    elif band == 'plw':
        b = 2

    #Get tiling info
    with open(output_folder+Master_filename, "rb") as f:
        obj = pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']

    #Make map array
    pred_map = np.empty_like(im)
    pred_map[:,:] = 0.0

    for i in range(start_tile, len(tiles)):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)-1))

        #Open tile
        tile = tiles[i]
        try:
            fi = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', "rb")
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                til = pickle.load(fi)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                fi.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
            else:
                fi.close()
                tilepri = til[band]
                tilepos = til['posterior']

                #Reshape posterior
                (samp,chains,params) = tilepos.stan_fit.shape
                flattened_post = tilepos.stan_fit.reshape(samp*chains,params)
                fvec = flattened_post[2999,0:tilepri.nsrc+1]

                #Map stuff (not quite sure what it is *shrug*)
                f = coo_matrix((fvec, (range(0,tilepri.nsrc+1),np.zeros(tilepri.nsrc+1))), shape=(tilepri.nsrc+1, 1))
                A = coo_matrix((tilepri.amat_data, (tilepri.amat_row, tilepri.amat_col)), shape=(tilepri.snpix, tilepri.nsrc+1))
                rmap_temp = (A*f)

                #Find pixels in tile
                ra, dec = w_pri.wcs_pix2world(tilepri.sx_pix, tilepri.sy_pix, 0)

                keep_pix = moc.sources_in_tile([tile], order, ra, dec)
                keep_pix = np.array(keep_pix)

                sx = tilepri.sx_pix * keep_pix
                sy = tilepri.sy_pix * keep_pix

                #Add tile to map
                pred_map[sy,sx] = np.asarray(rmap_temp.todense()).reshape(-1) +\
                                                          np.random.randn(tilepri.snpix)*\
                                                          np.sqrt(np.square(tilepri.snim)) #+np.square(conf_noise))
                #Reset variables in a vain attempt to reduce memory load
                tile = []
                til = []
                tilepri = []
                tilepos = []
                ra = []
                dec = []
                keep_pix = []
                sx = []
                sy = []

    fits.writeto(output_folder+'Master_Map_'+band+'.fits', pred_map, imhdu, overwrite = True)

def make_replicated_map_HEALpix(output_folder, Master_filename, tile_file_name, band, im, imhdu, w_pri, start_tile=0, noise=False):
    """Create a replicated map from XID+ with HEALpix tiles
       Will overwrite existing Master_Map_###.fits if one exists"""

    import pickle
    from scipy.sparse import coo_matrix
    from astropy import wcs
    from xidplus import moc_routines as moc
    from xidplus import posterior_maps as postmaps

    sam = 2999

    if band == 'psw':
        b = 0
    elif band == 'pmw':
        b = 1
    elif band == 'plw':
        b = 2
    elif band == 'mips24':
        b = 0
    elif band == 'green':
        b = 0
    elif band == 'red':
        b = 1
    else:
        print('Incompatible band entered')
        return

    #Get tiling info
    with open(output_folder+Master_filename, "rb") as f:
        obj = pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']

    #Make map array
    pred_map = np.empty_like(im)
    pred_map[:,:] = 0.0

    for i in range(start_tile, len(tiles)):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)-1))

        #Open tile
        tile = tiles[i]
        try:
            fi = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', "rb")
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                til = pickle.load(fi)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                fi.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
            else:
                fi.close()
                tilepri = til[band]
                tilepos = til['posterior']

                #Reshape posterior
                #(samp,chains,params) = tilepos.stan_fit.shape
                #flattened_post = tilepos.stan_fit.reshape(samp*chains,params)
                #fvec = flattened_post[2999,0:tilepri.nsrc+1]
                #
                #Map stuff (not quite sure what it is *shrug*)
                #f = coo_matrix((fvec, (range(0,tilepri.nsrc+1),np.zeros(tilepri.nsrc+1))), shape=(tilepri.nsrc+1, 1))
                #A = coo_matrix((tilepri.amat_data, (tilepri.amat_row, tilepri.amat_col)), shape=(tilepri.snpix, tilepri.nsrc+1))
                #rmap_temp = (A*f)

                rmap_array = postmaps.ymod_map(tilepri, tilepos.samples['src_f'][sam,b,:]).reshape(-1)

                #Find pixels in tile
                ra, dec = w_pri.wcs_pix2world(tilepri.sx_pix, tilepri.sy_pix, 0)

                keep_pix = moc.sources_in_tile([tile], order, ra, dec)
                keep_pix = np.array(keep_pix)

                sx = tilepri.sx_pix * keep_pix
                sy = tilepri.sy_pix * keep_pix

                #med_flux=tilepos.quantileGet(50)
                #bkg = np.full(tilepos.nsrc,med_flux[((1+mod)*tilepos.nsrc)+mod])[0]  #all background values are the same so just pick one
                #conf_noise = np.full(tilepos.nsrc,med_flux[(mod*tilepos.nsrc)+3+mod])[0] #all conf_noise is the same so just pick one
                if band == 'mips24':
                    bkg = tilepos.samples['bkg'][sam]
                else:
                    bkg = tilepos.samples['bkg'][sam,b]
                #Add tile to map
                if noise:
                    total_noise = np.random.normal(scale=np.sqrt(priors[b].snim**2+tilepos.samples['sigma_conf'][sam,b]**2))
                    pred_map[sy,sx] = rmap_array + bkg + total_noise
                    #                                          np.random.randn(tilepri.snpix)*\
                    #                                          np.sqrt(np.square(tilepri.snim)+np.square(conf_noise)) +\
                    #                                          bkg
                else:
                    #print(np.asarray(rmap_temp.todense()).reshape(-1))
                    #print(bkg)
                    #print(np.asarray(rmap_temp.todense()).reshape(-1) + bkg)
                    #pred_map[sy,sx] = np.asarray(rmap_temp.todense()).reshape(-1) + bkg
                    pred_map[sy,sx] = rmap_array + bkg

                #Reset variables in a vain attempt to reduce memory load
                tile = []
                til = []
                tilepri = []
                tilepos = []
                ra = []
                dec = []
                keep_pix = []
                sx = []
                sy = []

    fits.writeto(output_folder+'Replicated_Map_'+band+'.fits', pred_map, imhdu, overwrite = True)


def create_XIDp_MIPScat_from_Table(dataTable, prior24):
    import datetime
    from astropy.table import Table

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='ID', format='100A', array=dataTable['HELP_ID'])
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=dataTable['RA'])
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=dataTable['Dec'])
    c4 = fits.Column(name='F_MIPS_24', format='E', unit='uJy', array=dataTable['F_MIPS_24'])
    c5 = fits.Column(name='FErr_MIPS_24_u', format='E', unit='uJy', array=dataTable['FErr_MIPS_24_u'])
    c6 = fits.Column(name='FErr_MIPS_24_l', format='E', unit='uJy', array=dataTable['FErr_MIPS_24_l'])
    c7 = fits.Column(name='Bkg_MIPS_24', format='E', unit='MJy/sr', array=dataTable['Bkg_MIPS_24'])
    c8 = fits.Column(name='Sig_conf_MIPS_24', format='E',unit='MJy/sr', array=dataTable['Sig_conf_MIPS_24'])
    c9 = fits.Column(name='Rhat_MIPS_24', format='E', array=dataTable['Rhat_MIPS_24'])
    c10 = fits.Column(name='n_eff_MIPS_24', format='E', array=dataTable['n_eff_MIPS_24'])
    c11 = fits.Column(name='Pval_res_24', format='E', array=dataTable['Pval_res_24'])
    c12 = fits.Column(name='tile', format='K', array=dataTable['tile'].astype(float).astype(int))

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])

    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2')

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3')

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')
    tbhdu.header.set('TDESC4','24 Flux (at 50th percentile)',after='TUCD4')

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')
    tbhdu.header.set('TDESC5','24 Flux (at 84.1 percentile) ',after='TUCD5')

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')
    tbhdu.header.set('TDESC6','24 Flux (at 15.9 percentile)',after='TUCD6')

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')
    tbhdu.header.set('TDESC7','24 background',after='TUCD7')

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')
    tbhdu.header.set('TDESC8','24 residual confusion noise',after='TUCD8')

    tbhdu.header.set('TUCD9','stat.value',after='TFORM9')
    tbhdu.header.set('TDESC9','24 MCMC Convergence statistic',after='TUCD9')

    tbhdu.header.set('TUCD10','stat.value',after='TFORM10')
    tbhdu.header.set('TDESC10','24 MCMC independence statistic',after='TUCD10')
    
    tbhdu.header.set('TUCD11','stat.value',after='TFORM11')
    tbhdu.header.set('TDESC11','24 Bayes Pval residual statistic',after='TUCD11')
    
    tbhdu.header.set('TUCD12','value',after='TFORM12')
    tbhdu.header.set('TDESC12','Tile Number',after='TUCD12')

    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Pri_Cat'] = prior24.prior_cat_file
    prihdr['TITLE']   = 'MIPS XID+ catalogue'
    #prihdr['OBJECT']  = prior250.imphdu['OBJECT'] #I need to think if this needs to change
    prihdr['CREATOR'] = 'WP5'
    #prihdr['XID+V'] = git_version()
    prihdr['DATE']    = datetime.datetime.now().isoformat()
    prihdu = fits.PrimaryHDU(header=prihdr)

    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist

def make_master_MIPS_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=0):
    """Make master source catalogue from XID+ with HEALpix
       Master_filename needs the .pkl extension
       Creates Master_Catalogue.fits in output_folder
       Will overwrite existing Master_Catalogue.fits if one exists"""

    import pickle
    import dill
    from astropy.table import Table, vstack
    from xidplus import moc_routines, catalogue

    with open(output_folder+Master_filename, "rb") as f:
        obj=pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']
    prior24 = obj['psw']

    #Need a master catalogue to add to so do the 0th tile outside the loop
    for i in range(start_tile, len(tiles)):
        print('On tile '+str(i)+' of ' + str(len(tiles)))

        tile = tiles[i]
        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing?')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data?')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()
                main_start = i+1
                break

    tilemips = obj['psw']
    tilepos = obj['posterior']

    master_cat = catalogue.create_MIPS_cat(tilepos, tilemips)
    #Only keep sources inside the tile
    kept_sources = moc_routines.sources_in_tile([tile], order, tilemips.sra, tilemips.sdec)
    kept_sources = np.array(kept_sources)

    master_cat[1].data = master_cat[1].data[kept_sources]
    #Create the master catalogue table
    dataTable = Table(master_cat[1].data)
    dataTable['tile'] = main_start-1

    for i in range(main_start, len(tiles)):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)))
        tile = tiles[i]

        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()

                tilemips = obj['psw']
                tilepos = obj['posterior']

                hdulist = catalogue.create_MIPS_cat(tilepos, tilemips)
                #Only keep sources inside the tile
                kept_sources = moc_routines.sources_in_tile([tile], order, tilemips.sra, tilemips.sdec)
                kept_sources = np.array(kept_sources)

                hdulist[1].data = hdulist[1].data[kept_sources]
                tileTable = Table(hdulist[1].data)
                tileTable['tile'] = np.full(len(tileTable),i)
                #Append the tile catalogue to the master catalogue
                dataTable = vstack([dataTable, tileTable])
                print(len(dataTable))

    #Write the master catalogue to fits
    master_cat = create_XIDp_MIPScat_from_Table(dataTable, prior24)
    master_cat.writeto(output_folder+'Master_Catalogue_MIPS24.fits', overwrite = True)
    
def create_MIPS_cat_post(posterior,prior24):
    """
    Create SPIRE catalogue of posterior

    :param posterior: MIPS xidplus.posterior class
    :param prior24: MIPS 24 xidplus.prior class
    :return: fits hdulist
    """
    import datetime
    nsrc=posterior.nsrc
    samp_chains=len(posterior.samples['src_f'])

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='HELP_ID', format='100A', array=prior24.ID)
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=prior24.sra)
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=prior24.sdec)
    c4 = fits.Column(name='F_MIPS_24', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],50.0,axis=0))
    c5 = fits.Column(name='FErr_MIPS_24_u', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],84.1,axis=0))
    c6 = fits.Column(name='FErr_MIPS_24_l', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],15.9,axis=0))
    c7 = fits.Column(name='Bkg_MIPS_24', format='E', unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['bkg'][:,0],50.0,axis=0)))
    c8 = fits.Column(name='Sig_conf_MIPS_24', format='E',unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['sigma_conf'][:,0],50.0,axis=0)))
    c9 = fits.Column(name='Rhat_MIPS_24', format='E', array=posterior.Rhat['src_f'][:,0])
    c10 = fits.Column(name='n_eff_MIPS_24', format='E', array=posterior.n_eff['src_f'][:,0])
    c11 = fits.Column(name='Post_MIPS_24', format=str(samp_chains)+'E', array=posterior.samples['src_f'][:,0,:].T)
    #c12 = fits.Column(name='Post_bkg_MIPS_24, format=str(samp_chains)+'E', unit='mJy', array=np.repeat(posterior.samples['bkg'][:,0], nsrc, axis=0))

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11]) #,c12])
    
    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2')

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3')

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')
    tbhdu.header.set('TDESC4','24 Flux (at 50th percentile)',after='TUCD4')

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')
    tbhdu.header.set('TDESC5','24 Flux (at 84.1 percentile) ',after='TUCD5')

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')
    tbhdu.header.set('TDESC6','24 Flux (at 15.9 percentile)',after='TUCD6')

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')
    tbhdu.header.set('TDESC7','24 background',after='TUCD7')

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')
    tbhdu.header.set('TDESC8','24 residual confusion noise',after='TUCD8')

    tbhdu.header.set('TUCD9','stat.value',after='TFORM9')
    tbhdu.header.set('TDESC9','24 MCMC Convergence statistic',after='TUCD9')

    tbhdu.header.set('TUCD10','stat.value',after='TFORM10')
    tbhdu.header.set('TDESC10','24 MCMC independence statistic',after='TUCD10')

    tbhdu.header.set('TUCD11','phot.flux.density',after='TFORM11')
    tbhdu.header.set('TDESC11','24 samples',after='TUCD11')

    #tbhdu.header.set('TUCD12','phot.flux.density',after='TFORM12')
    #tbhdu.header.set('TDESC12','24 bkg samples',after='TUCD12')


    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Prior_C'] = prior24.prior_cat_file
    prihdr['TITLE']   = 'MIPS XID catalogue'        
    #prihdr['OBJECT']  = prior24.imphdu['OBJECT'] #I need to think if this needs to change                              
    prihdr['CREATOR'] = 'WP5'                                 
    #prihdr['VERSION'] = 'beta'                                 
    prihdr['DATE']    = datetime.datetime.now().isoformat()              
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist
    
    
def create_MIPS_cat_post_from_Table(dataTable, prior24, samples_chains):
    """creates the XIDp catalogue in fits format required by HeDaM"""
    import datetime
    from astropy.table import Table

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='ID', format='100A', array=dataTable['HELP_ID'])
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=dataTable['RA'])
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=dataTable['Dec'])
    c4 = fits.Column(name='F_MIPS_24', format='E', unit='uJy', array=dataTable['F_MIPS_24'])
    c5 = fits.Column(name='FErr_MIPS_24_u', format='E', unit='uJy', array=dataTable['FErr_MIPS_24_u'])
    c6 = fits.Column(name='FErr_MIPS_24_l', format='E', unit='uJy', array=dataTable['FErr_MIPS_24_l'])
    c7 = fits.Column(name='Bkg_MIPS_24', format='E', unit='uJy', array=dataTable['Bkg_MIPS_24'])
    c8 = fits.Column(name='Sig_conf_MIPS_24', format='E',unit='MJy/Beam', array=dataTable['Sig_conf_MIPS_24'])
    c9 = fits.Column(name='Rhat_MIPS_24', format='E', array=dataTable['Rhat_MIPS_24'])
    c10 = fits.Column(name='n_eff_MIPS_24', format='E', array=dataTable['n_eff_MIPS_24'])
    c11 = fits.Column(name='Post_MIPS_24', format=str(samples_chains)+'E', unit='uJy', array=dataTable['Post_MIPS_24'])
    #c12 = fits.Column(name='Post_bkg_MIPS_24', format=str(samples*chains)+'E', unit='mJy', array=dataTable['Post_bkg_MIPS_24])

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11]) #,c12])

    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2')

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3')

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')
    tbhdu.header.set('TDESC4','24 Flux (at 50th percentile)',after='TUCD4')

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')
    tbhdu.header.set('TDESC5','24 Flux (at 84.1 percentile) ',after='TUCD5')

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')
    tbhdu.header.set('TDESC6','24 Flux (at 15.9 percentile)',after='TUCD6')

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')
    tbhdu.header.set('TDESC7','24 background',after='TUCD7')

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')
    tbhdu.header.set('TDESC8','24 residual confusion noise',after='TUCD8')

    tbhdu.header.set('TUCD9','stat.value',after='TFORM9')
    tbhdu.header.set('TDESC9','24 MCMC Convergence statistic',after='TUCD9')

    tbhdu.header.set('TUCD10','stat.value',after='TFORM10')
    tbhdu.header.set('TDESC10','24 MCMC independence statistic',after='TUCD10')

    tbhdu.header.set('TUCD11','phot.flux.density',after='TFORM11')
    tbhdu.header.set('TDESC11','24 samples',after='TUCD11')

    #tbhdu.header.set('TUCD12','phot.flux.density',after='TFORM12')
    #tbhdu.header.set('TDESC12','24 bkg samples',after='TUCD12')
    

    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Pri_Cat'] = prior24.prior_cat_file
    prihdr['TITLE']   = 'MIPS XID+ Posterior catalogue'
    #prihdr['OBJECT']  = prior24.imphdu['OBJECT'] #I need to think if this needs to change
    prihdr['CREATOR'] = 'WP5'
    #prihdr['XID+V'] = git_version()
    prihdr['DATE']    = datetime.datetime.now().isoformat()
    prihdu = fits.PrimaryHDU(header=prihdr)

    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist

def make_master_MIPS_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=0, end_tile=0):
    """Make master source catalogue from XID+ with HEALpix
       Master_filename needs the .pkl extension
       Creates Master_Catalogue.fits in output_folder
       Will overwrite existing Master_Catalogue.fits if one exists"""

    import pickle
    import dill
    from astropy.table import Table, vstack
    from xidplus import moc_routines, catalogue

    with open(output_folder+Master_filename, "rb") as f:
        obj=pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']
    prior24 = obj['psw']

    if end_tile == 0:
        end_tile = len(tiles)

    #Need a master catalogue to add to so do the 0th tile outside the loop
    for i in range(start_tile, end_tile):
        print('On tile '+str(i)+' of ' + str(len(tiles)))

        tile = tiles[i]
        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing?')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data?')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()
                main_start = i+1
                break

    tilepsw = obj['psw']
    tilepos = obj['posterior']

    master_cat = create_MIPS_cat_post(tilepos, tilepsw)
    #Only keep sources inside the tile
    kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
    kept_sources = np.array(kept_sources)

    master_cat[1].data = master_cat[1].data[kept_sources]
    #Create the master catalogue table
    dataTable = Table(master_cat[1].data)
    dataTable['tile'] = main_start-1

    for i in range(main_start, end_tile):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)))
        tile = tiles[i]

        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()

                tilepsw = obj['psw']
                tilepos = obj['posterior']

                hdulist = create_MIPS_cat_post(tilepos, tilepsw)
                #Only keep sources inside the tile
                kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
                kept_sources = np.array(kept_sources)

                hdulist[1].data = hdulist[1].data[kept_sources]
                tileTable = Table(hdulist[1].data)
                #Append the tile catalogue to the master catalogue
                dataTable = vstack([dataTable, tileTable])

    #Write the master catalogue to fits
    samples_chains = len(tilepos.samples['src_f'])
    master_cat = create_MIPS_cat_post_from_Table(dataTable, prior24, samples_chains)
    master_cat.writeto(output_folder+'Master_Post_Catalogue_MIPS_'+str(start_tile)+'.fits', overwrite = True)
    

def create_XIDp_PACScat_from_Table(dataTable, prior100):
    """creates the XIDp catalogue in fits format required by HeDaM"""
    import datetime
    from astropy.table import Table

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='ID', format='100A', array=dataTable['HELP_ID'])
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=dataTable['RA'])
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=dataTable['Dec'])
    c4 = fits.Column(name='F_PACS_100', format='E', unit='mJy', array=dataTable['F_PACS_100'])
    c5 = fits.Column(name='FErr_PACS_100_u', format='E', unit='mJy', array=dataTable['FErr_PACS_100_u'])
    c6 = fits.Column(name='FErr_PACS_100_l', format='E', unit='mJy', array=dataTable['FErr_PACS_100_l'])
    c7 = fits.Column(name='F_PACS_160', format='E', unit='mJy', array=dataTable['F_PACS_160'])
    c8 = fits.Column(name='FErr_PACS_160_u', format='E', unit='mJy', array=dataTable['FErr_PACS_160_u'])
    c9 = fits.Column(name='FErr_PACS_160_l', format='E', unit='mJy', array=dataTable['FErr_PACS_160_l'])
    c10 = fits.Column(name='Bkg_PACS_100', format='E', unit='mJy', array=dataTable['Bkg_PACS_100'])
    c11 = fits.Column(name='Bkg_PACS_160', format='E', unit='mJy', array=dataTable['Bkg_PACS_160'])
    c12 = fits.Column(name='Sig_conf_PACS_100', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_PACS_100'])
    c13 = fits.Column(name='Sig_conf_PACS_160', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_PACS_160'])
    c14 = fits.Column(name='Rhat_PACS_100', format='E', array=dataTable['Rhat_PACS_100'])
    c15 = fits.Column(name='Rhat_PACS_160', format='E', array=dataTable['Rhat_PACS_160'])
    c16 = fits.Column(name='n_eff_PACS_100', format='E', array=dataTable['n_eff_PACS_100'])
    c17 = fits.Column(name='n_eff_PACS_160', format='E', array=dataTable['n_eff_PACS_160'])
    c18 = fits.Column(name='Pval_res_100', format='E', array=dataTable['Pval_res_100'])
    c19 = fits.Column(name='Pval_res_160', format='E', array=dataTable['Pval_res_100'])
    c20 = fits.Column(name='tile', format='K', array=dataTable['tile'].astype(float).astype(int))

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20])

    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2')

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3')

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')
    tbhdu.header.set('TDESC4','100 Flux (at 50th percentile)',after='TUCD4')

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')
    tbhdu.header.set('TDESC5','100 Flux (at 84.1 percentile) ',after='TUCD5')

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')
    tbhdu.header.set('TDESC6','100 Flux (at 15.9 percentile)',after='TUCD6')

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')
    tbhdu.header.set('TDESC7','160 Flux (at 50th percentile)',after='TUCD7')

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')
    tbhdu.header.set('TDESC8','160 Flux (at 84.1 percentile) ',after='TUCD8')

    tbhdu.header.set('TUCD9','phot.flux.density',after='TUNIT9')
    tbhdu.header.set('TDESC9','160 Flux (at 15.9 percentile)',after='TUCD9')

    tbhdu.header.set('TUCD10','phot.flux.density',after='TUNIT10')
    tbhdu.header.set('TDESC10','100 background',after='TUCD10')

    tbhdu.header.set('TUCD11','phot.flux.density',after='TUNIT11')
    tbhdu.header.set('TDESC11','160 background',after='TUCD11')

    tbhdu.header.set('TUCD12','phot.flux.density',after='TUNIT12')
    tbhdu.header.set('TDESC12','100 residual confusion noise',after='TUCD12')

    tbhdu.header.set('TUCD13','phot.flux.density',after='TUNIT13')
    tbhdu.header.set('TDESC13','160 residual confusion noise',after='TUCD13')

    tbhdu.header.set('TUCD14','stat.value',after='TFORM14')
    tbhdu.header.set('TDESC14','100 MCMC Convergence statistic',after='TUCD14')

    tbhdu.header.set('TUCD15','stat.value',after='TFORM15')
    tbhdu.header.set('TDESC15','160 MCMC Convergence statistic',after='TUCD15')

    tbhdu.header.set('TUCD16','stat.value',after='TFORM16')
    tbhdu.header.set('TDESC16','100 MCMC independence statistic',after='TUCD16')

    tbhdu.header.set('TUCD17','stat.value',after='TFORM17')
    tbhdu.header.set('TDESC17','160 MCMC independence statistic',after='TUCD17')
    
    tbhdu.header.set('TUCD18','stat.value',after='TFORM18')
    tbhdu.header.set('TDESC18','100 Bayes Pval residual statistic',after='TUCD18')

    tbhdu.header.set('TUCD19','stat.value',after='TFORM19')
    tbhdu.header.set('TDESC19','160 Bayes Pval residual statistic',after='TUCD19')
    
    tbhdu.header.set('TUCD20','value',after='TFORM20')
    tbhdu.header.set('TDESC20','Tile Number',after='TUCD20')

    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Pri_Cat'] = prior100.prior_cat_file
    prihdr['TITLE']   = 'PACS XID+ catalogue'
    #prihdr['OBJECT']  = prior250.imphdu['OBJECT'] #I need to think if this needs to change
    prihdr['CREATOR'] = 'WP5'
    #prihdr['XID+V'] = git_version()
    prihdr['DATE']    = datetime.datetime.now().isoformat()
    prihdu = fits.PrimaryHDU(header=prihdr)

    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist

def make_master_PACS_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=0):
    """Make master source catalogue from XID+ with HEALpix
       Master_filename needs the .pkl extension
       Creates Master_Catalogue.fits in output_folder
       Will overwrite existing Master_Catalogue.fits if one exists"""

    import pickle
    import dill
    from astropy.table import Table, vstack
    from xidplus import moc_routines, catalogue

    with open(output_folder+Master_filename, "rb") as f:
        obj=pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']
    prior100 = obj['green']

    #Need a master catalogue to add to so do the 0th tile outside the loop
    for i in range(start_tile, len(tiles)):
        print('On tile '+str(i)+' of ' + str(len(tiles)))

        tile = tiles[i]
        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing?')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data?')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()
                main_start = i+1
                break

    tilegrn = obj['green']
    tilered = obj['red']
    tilepos = obj['posterior']

    master_cat = catalogue.create_PACS_cat(tilepos, tilegrn, tilered)
    #Only keep sources inside the tile
    kept_sources = moc_routines.sources_in_tile([tile], order, tilegrn.sra, tilegrn.sdec)
    kept_sources = np.array(kept_sources)

    master_cat[1].data = master_cat[1].data[kept_sources]
    #Create the master catalogue table
    dataTable = Table(master_cat[1].data)
    dataTable['tile'] = main_start-1

    for i in range(main_start, len(tiles)):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)))
        tile = tiles[i]

        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()

                tilegrn = obj['green']
                tilered = obj['red']
                tilepos = obj['posterior']

                hdulist = catalogue.create_PACS_cat(tilepos, tilegrn, tilered)
                #Only keep sources inside the tile
                kept_sources = moc_routines.sources_in_tile([tile], order, tilegrn.sra, tilegrn.sdec)
                kept_sources = np.array(kept_sources)

                hdulist[1].data = hdulist[1].data[kept_sources]
                tileTable = Table(hdulist[1].data)
                tileTable['tile'] = np.full(len(tileTable),i)
                #Append the tile catalogue to the master catalogue
                dataTable = vstack([dataTable, tileTable])

    #Write the master catalogue to fits
    master_cat = create_XIDp_PACScat_from_Table(dataTable, prior100)
    master_cat.writeto(output_folder+'Master_Catalogue_PACS.fits', overwrite = True)

def create_PACS_cat_post(posterior,prior100,prior160):
    """
    Create SPIRE catalogue of posterior

    :param posterior: PACS xidplus.posterior class
    :param prior100: PACS 100 xidplus.prior class
    :param prior160: PACS 160 xidplus.prior class
    :return: fits hdulist
    """
    import datetime
    nsrc=posterior.nsrc
    samp_chains=len(posterior.samples['src_f'])

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='HELP_ID', format='100A', array=prior100.ID)
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=prior100.sra)
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=prior100.sdec)
    c4 = fits.Column(name='F_PACS_100', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],50.0,axis=0))
    c5 = fits.Column(name='FErr_PACS_100_u', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],84.1,axis=0))
    c6 = fits.Column(name='FErr_PACS_100_l', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,0,:],15.9,axis=0))
    c7 = fits.Column(name='F_PACS_160', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,1,:],50.0,axis=0))
    c8 = fits.Column(name='FErr_PACS_160_u', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,1,:],84.1,axis=0))
    c9 = fits.Column(name='FErr_PACS_160_l', format='E', unit='mJy',
                     array=np.percentile(posterior.samples['src_f'][:,1,:],15.9,axis=0))
    c10 = fits.Column(name='Bkg_PACS_100', format='E', unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['bkg'][:,0],50.0,axis=0)))
    c11 = fits.Column(name='Bkg_PACS_160', format='E', unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['bkg'][:,1],50.0,axis=0)))
    c12 = fits.Column(name='Sig_conf_PACS_100', format='E',unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['sigma_conf'][:,0],50.0,axis=0)))
    c13 = fits.Column(name='Sig_conf_PACS_160', format='E',unit='mJy/Beam',
                      array=np.full(nsrc,np.percentile(posterior.samples['sigma_conf'][:,1],50.0,axis=0)))
    c14 = fits.Column(name='Rhat_PACS_100', format='E', array=posterior.Rhat['src_f'][:,0])
    c15 = fits.Column(name='Rhat_PACS_160', format='E', array=posterior.Rhat['src_f'][:,1])
    c16 = fits.Column(name='n_eff_PACS_100', format='E', array=posterior.n_eff['src_f'][:,0])
    c17 = fits.Column(name='n_eff_PACS_160', format='E', array=posterior.n_eff['src_f'][:,1])
    c18 = fits.Column(name='Post_PACS_100', format=str(samp_chains)+'E', array=posterior.samples['src_f'][:,0,:].T)
    c19 = fits.Column(name='Post_PACS_160', format=str(samp_chains)+'E', array=posterior.samples['src_f'][:,1,:].T)
    #c20 = fits.Column(name='Post_bkg_PACS_100', format=str(samp_chains)+'E', unit='mJy', array=np.repeat(posterior.samples['bkg'][:,0], nsrc, axis=0))
    #c21 = fits.Column(name='Post_bkg_PACS_160', format=str(samp_chains)+'E', unit='mJy', array=np.repeat(posterior.samples['bkg'][:,1], nsrc, axis=0))

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20]) #,c21])
    
    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2')

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3')

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')
    tbhdu.header.set('TDESC4','100 Flux (at 50th percentile)',after='TUCD4')

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')
    tbhdu.header.set('TDESC5','100 Flux (at 84.1 percentile) ',after='TUCD5')

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')
    tbhdu.header.set('TDESC6','100 Flux (at 15.9 percentile)',after='TUCD6')

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')
    tbhdu.header.set('TDESC7','160 Flux (at 50th percentile)',after='TUCD7')

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')
    tbhdu.header.set('TDESC8','160 Flux (at 84.1 percentile) ',after='TUCD8')

    tbhdu.header.set('TUCD9','phot.flux.density',after='TUNIT9')
    tbhdu.header.set('TDESC9','160 Flux (at 15.9 percentile)',after='TUCD9')

    tbhdu.header.set('TUCD10','phot.flux.density',after='TUNIT10')
    tbhdu.header.set('TDESC10','100 background',after='TUCD10')

    tbhdu.header.set('TUCD11','phot.flux.density',after='TUNIT11')
    tbhdu.header.set('TDESC11','160 background',after='TUCD11')

    tbhdu.header.set('TUCD12','phot.flux.density',after='TUNIT12')
    tbhdu.header.set('TDESC12','100 residual confusion noise',after='TUCD12')

    tbhdu.header.set('TUCD13','phot.flux.density',after='TUNIT13')
    tbhdu.header.set('TDESC13','160 residual confusion noise',after='TUCD13')

    tbhdu.header.set('TUCD14','stat.value',after='TFORM14')
    tbhdu.header.set('TDESC14','100 MCMC Convergence statistic',after='TUCD14')

    tbhdu.header.set('TUCD15','stat.value',after='TFORM15')
    tbhdu.header.set('TDESC15','160 MCMC Convergence statistic',after='TUCD15')

    tbhdu.header.set('TUCD16','stat.value',after='TFORM16')
    tbhdu.header.set('TDESC16','100 MCMC independence statistic',after='TUCD16')

    tbhdu.header.set('TUCD17','stat.value',after='TFORM17')
    tbhdu.header.set('TDESC17','160 MCMC independence statistic',after='TUCD17')

    tbhdu.header.set('TUCD18','phot.flux.density',after='TFORM18')
    tbhdu.header.set('TDESC18','100 samples',after='TUCD18')

    tbhdu.header.set('TUCD19','phot.flux.density',after='TFORM19')
    tbhdu.header.set('TDESC19','160 samples',after='TUCD19')

    #tbhdu.header.set('TUCD20','phot.flux.density',after='TFORM20')
    #tbhdu.header.set('TDESC20','100 bkg samples',after='TUCD20')

    #tbhdu.header.set('TUCD21','phot.flux.density',after='TFORM21')
    #tbhdu.header.set('TDESC21','160 bkg samples',after='TUCD21')


    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Prior_C'] = prior100.prior_cat_file
    prihdr['TITLE']   = 'PACS XID catalogue'        
    #prihdr['OBJECT']  = prior100.imphdu['OBJECT'] #I need to think if this needs to change                              
    prihdr['CREATOR'] = 'WP5'                                 
    #prihdr['VERSION'] = 'beta'                                 
    prihdr['DATE']    = datetime.datetime.now().isoformat()              
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist
    
def create_PACS_cat_post_from_Table(dataTable, prior100, samples_chains):
    """creates the XIDp catalogue in fits format required by HeDaM"""
    import datetime
    from astropy.table import Table

    #----table info-----------------------
    #first define columns
    c1 = fits.Column(name='ID', format='100A', array=dataTable['HELP_ID'])
    c2 = fits.Column(name='RA', format='D', unit='degrees', array=dataTable['RA'])
    c3 = fits.Column(name='Dec', format='D', unit='degrees', array=dataTable['Dec'])
    c4 = fits.Column(name='F_PACS_100', format='E', unit='mJy', array=dataTable['F_PACS_100'])
    c5 = fits.Column(name='FErr_PACS_100_u', format='E', unit='mJy', array=dataTable['FErr_PACS_100_u'])
    c6 = fits.Column(name='FErr_PACS_100_l', format='E', unit='mJy', array=dataTable['FErr_PACS_100_l'])
    c7 = fits.Column(name='F_PACS_160', format='E', unit='mJy', array=dataTable['F_PACS_160'])
    c8 = fits.Column(name='FErr_PACS_160_u', format='E', unit='mJy', array=dataTable['FErr_PACS_160_u'])
    c9 = fits.Column(name='FErr_PACS_160_l', format='E', unit='mJy', array=dataTable['FErr_PACS_160_l'])
    c10 = fits.Column(name='Bkg_PACS_100', format='E', unit='mJy', array=dataTable['Bkg_PACS_100'])
    c11 = fits.Column(name='Bkg_PACS_160', format='E', unit='mJy', array=dataTable['Bkg_PACS_160'])
    c12 = fits.Column(name='Sig_conf_PACS_100', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_PACS_100'])
    c13 = fits.Column(name='Sig_conf_PACS_160', format='E',unit='mJy/Beam', array=dataTable['Sig_conf_PACS_160'])
    c14 = fits.Column(name='Rhat_PACS_100', format='E', array=dataTable['Rhat_PACS_100'])
    c15 = fits.Column(name='Rhat_PACS_160', format='E', array=dataTable['Rhat_PACS_160'])
    c16 = fits.Column(name='n_eff_PACS_100', format='E', array=dataTable['n_eff_PACS_100'])
    c17 = fits.Column(name='n_eff_PACS_160', format='E', array=dataTable['n_eff_PACS_160'])
    c18 = fits.Column(name='Post_PACS_100', format=str(samples_chains)+'E', unit='mJy', array=dataTable['Post_PACS_100'])
    c19 = fits.Column(name='Post_PACS_160', format=str(samples_chains)+'E', unit='mJy', array=dataTable['Post_PACS_160'])
    #c20 = fits.Column(name='Post_bkg_PACS_100', format=str(samples*chains)+'E', unit='mJy', array=dataTable['Post_bkg_PACS_100'])
    #c21 = fits.Column(name='Post_bkg_PACS_100', format=str(samples*chains)+'E', unit='mJy', array=dataTable['Post_bkg_PACS_100'])

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19]) #,c20,c21])

    tbhdu.header.set('TUCD1','ID',after='TFORM1')
    tbhdu.header.set('TDESC1','ID of source',after='TUCD1')

    tbhdu.header.set('TUCD2','pos.eq.RA',after='TUNIT2')
    tbhdu.header.set('TDESC2','R.A. of object J2000',after='TUCD2')

    tbhdu.header.set('TUCD3','pos.eq.DEC',after='TUNIT3')
    tbhdu.header.set('TDESC3','Dec. of object J2000',after='TUCD3')

    tbhdu.header.set('TUCD4','phot.flux.density',after='TUNIT4')
    tbhdu.header.set('TDESC4','100 Flux (at 50th percentile)',after='TUCD4')

    tbhdu.header.set('TUCD5','phot.flux.density',after='TUNIT5')
    tbhdu.header.set('TDESC5','100 Flux (at 84.1 percentile) ',after='TUCD5')

    tbhdu.header.set('TUCD6','phot.flux.density',after='TUNIT6')
    tbhdu.header.set('TDESC6','100 Flux (at 15.9 percentile)',after='TUCD6')

    tbhdu.header.set('TUCD7','phot.flux.density',after='TUNIT7')
    tbhdu.header.set('TDESC7','160 Flux (at 50th percentile)',after='TUCD7')

    tbhdu.header.set('TUCD8','phot.flux.density',after='TUNIT8')
    tbhdu.header.set('TDESC8','160 Flux (at 84.1 percentile) ',after='TUCD8')

    tbhdu.header.set('TUCD9','phot.flux.density',after='TUNIT9')
    tbhdu.header.set('TDESC9','160 Flux (at 15.9 percentile)',after='TUCD9')

    tbhdu.header.set('TUCD10','phot.flux.density',after='TUNIT10')
    tbhdu.header.set('TDESC10','100 background',after='TUCD10')

    tbhdu.header.set('TUCD11','phot.flux.density',after='TUNIT11')
    tbhdu.header.set('TDESC11','160 background',after='TUCD11')

    tbhdu.header.set('TUCD12','phot.flux.density',after='TUNIT12')
    tbhdu.header.set('TDESC12','100 residual confusion noise',after='TUCD12')

    tbhdu.header.set('TUCD13','phot.flux.density',after='TUNIT13')
    tbhdu.header.set('TDESC13','160 residual confusion noise',after='TUCD13')

    tbhdu.header.set('TUCD14','stat.value',after='TFORM14')
    tbhdu.header.set('TDESC14','100 MCMC Convergence statistic',after='TUCD14')

    tbhdu.header.set('TUCD15','stat.value',after='TFORM15')
    tbhdu.header.set('TDESC15','160 MCMC Convergence statistic',after='TUCD15')

    tbhdu.header.set('TUCD16','stat.value',after='TFORM16')
    tbhdu.header.set('TDESC16','100 MCMC independence statistic',after='TUCD16')

    tbhdu.header.set('TUCD17','stat.value',after='TFORM17')
    tbhdu.header.set('TDESC17','160 MCMC independence statistic',after='TUCD17')

    tbhdu.header.set('TUCD18','phot.flux.density',after='TFORM18')
    tbhdu.header.set('TDESC19','100 samples',after='TUCD19')

    tbhdu.header.set('TUCD20','phot.flux.density',after='TFORM20')
    tbhdu.header.set('TDESC20','160 samples',after='TUCD20')

    #tbhdu.header.set('TUCD21','phot.flux.density',after='TFORM21')
    #tbhdu.header.set('TDESC21','100 bkg samples',after='TUCD21')

    #tbhdu.header.set('TUCD22','phot.flux.density',after='TFORM22')
    #tbhdu.header.set('TDESC22','160 bkg samples',after='TUCD22')
    

    #----Primary header-----------------------------------
    prihdr = fits.Header()
    prihdr['Pri_Cat'] = prior100.prior_cat_file
    prihdr['TITLE']   = 'PACS XID+ Posterior catalogue'
    #prihdr['OBJECT']  = prior100.imphdu['OBJECT'] #I need to think if this needs to change
    prihdr['CREATOR'] = 'WP5'
    #prihdr['XID+V'] = git_version()
    prihdr['DATE']    = datetime.datetime.now().isoformat()
    prihdu = fits.PrimaryHDU(header=prihdr)

    thdulist = fits.HDUList([prihdu, tbhdu])
    return thdulist

def make_master_PACS_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=0, end_tile=0):
    """Make master source catalogue from XID+ with HEALpix
       Master_filename needs the .pkl extension
       Creates Master_Catalogue.fits in output_folder
       Will overwrite existing Master_Catalogue.fits if one exists"""

    import pickle
    import dill
    from astropy.table import Table, vstack
    from xidplus import moc_routines, catalogue

    with open(output_folder+Master_filename, "rb") as f:
        obj=pickle.load(f)
        f.close()

    order = obj['order']
    tiles = obj['tiles']
    prior100 = obj['psw']

    if end_tile == 0:
        end_tile = len(tiles)

    #Need a master catalogue to add to so do the 0th tile outside the loop
    for i in range(start_tile, end_tile):
        print('On tile '+str(i)+' of ' + str(len(tiles)))

        tile = tiles[i]
        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing?')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data?')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()
                main_start = i+1
                break

    tilepsw = obj['psw']
    tilepmw = obj['pmw']
    tilepos = obj['posterior']

    master_cat = create_PACS_cat_post(tilepos, tilepsw, tilepmw)
    #Only keep sources inside the tile
    kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
    kept_sources = np.array(kept_sources)

    master_cat[1].data = master_cat[1].data[kept_sources]
    #Create the master catalogue table
    dataTable = Table(master_cat[1].data)
    dataTable['tile'] = main_start-1

    for i in range(main_start, end_tile):
        print('On tile ' + str(i) + ' of ' + str(len(tiles)))
        tile = tiles[i]

        try:
            f = open(output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl', 'rb')
        except IOError:
            print('IOError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl missing')

        else:
            try:
                obj = pickle.load(f)
            except EOFError:
                print('EOFError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl has no data')
                f.close()
            except ValueError:
                print('ValueError: '+output_folder+tile_file_name+str(tile)+'_'+str(order)+'.pkl')
                f.close()

            else:
                f.close()

                tilepsw = obj['psw']
                tilepmw = obj['pmw']
                tilepos = obj['posterior']

                hdulist = fcreate_PACS_cat_post(tilepos, tilepsw, tilepmw)
                #Only keep sources inside the tile
                kept_sources = moc_routines.sources_in_tile([tile], order, tilepsw.sra, tilepsw.sdec)
                kept_sources = np.array(kept_sources)

                hdulist[1].data = hdulist[1].data[kept_sources]
                tileTable = Table(hdulist[1].data)
                #Append the tile catalogue to the master catalogue
                dataTable = vstack([dataTable, tileTable])

    #Write the master catalogue to fits
    samples_chains = len(tilepos.samples['src_f'])
    master_cat = create_PACS_cat_post_from_Table(dataTable, prior100, samples_chains)
    master_cat.writeto(output_folder+'Master_Post_Catalogue_PACS_'+str(start_tile)+'.fits', overwrite = True)
