from celery_worker import celery
from utils import substructure_search, logger
import asyncio

@celery.task(name='src.tasks.substructure_search_task')
def substructure_search_task(substructure: str, molecules_list: list):

    try:
        matches = asyncio.run(substructure_search(molecules_list, substructure))
        return matches
    except Exception as e:
        logger.error(f"An error occurred during substructure search: {str(e)}")
        raise e