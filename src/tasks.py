from celery import shared_task
from redis_utils import set_cache
from utils import substructure_search, logger
import asyncio


@shared_task(name='src.tasks.substructure_search_task')
def substructure_search_task(substructure: str, molecules_list: list):
    try:
        matches = asyncio.run(substructure_search(molecules_list, substructure))

        cache_key = f"subsearch:{substructure}"
        set_cache(cache_key, matches, expiration=3600)

        return matches
    except Exception as e:
        logger.error(f"An error occurred during substructure search: {str(e)}")
        raise e
