from igvf_client import Configuration
from igvf_client import ApiClient
from igvf_client import IgvfApi


def get_igvf_auth_and_api(igvf_site: str = 'production'):
    """Set up IGVF data portal access and set up IGVF python client api.
    """
    if igvf_site == 'production':
        host = 'https://api.data.igvf.org'
    elif igvf_site == 'sandbox':
        host = 'https://api.sandbox.igvf.org'
    config = Configuration(
        access_key=os.environ["IGVF_API_KEY"],
        secret_access_key=os.environ["IGVF_SECRET_KEY"],
        host=host
    )
    client = ApiClient(config)
    return IgvfApi(client)
