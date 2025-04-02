# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:00:26 2022

@author: Dustin Garner
"""

from fafbseg import flywire
import flycode.neuprint_reading as neuread


flywire.set_default_dataset("production")
try:
    neuread.initialize_client()
except:
    print('Unable to initialize Neuprint.')


"""
Get a token from by logging in with your email account that is registered with
Flywire from https://global.daf-apis.com/auth/api/v1/user/token.
Then run this function:
    
token = ""
flywire.set_chunkedgraph_secret(token)
"""

