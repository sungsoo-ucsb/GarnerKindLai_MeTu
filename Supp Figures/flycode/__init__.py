<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:00:26 2022

@author: dusti
"""

from fafbseg import flywire


flywire.set_default_dataset("production")

#flywire.set_chunkedgraph_secret("440eb6c6e88bc08792e8745dbffedff6")


=======
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:00:26 2022

@author: Dustin Garner
"""

from fafbseg import flywire


flywire.set_default_dataset("production")


"""
Get a token from by logging in with your email account that is registered with
Flywire from https://global.daf-apis.com/auth/api/v1/user/token.
Then run this function:
    
token = ""
flywire.set_chunkedgraph_secret(token)
"""

>>>>>>> Stashed changes
